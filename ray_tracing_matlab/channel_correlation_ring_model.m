clear all; close all;
%% Main Parameters
write = "no"
maxNoRef= 5             ;
txpower= 10;  % in dbm
f = 2e9; % in Hz
lambda = physconst('lightspeed')/f;
noLoc = 5;

%% Map
mapFileName = "one_ring_model.osm";
viewer = siteviewer("Buildings",mapFileName ,"Basemap","topographic");

%% Tranceivers
rxNames = {'Bob', 'Eve'};
assert(length(rxNames)==2,"This code is programmed for two receivers only")
Txroute  = [...
            51.456651,-2.585808;...% This the route followed by the tx;
            51.456483,-2.586995
            ];    

        
rxLocations =[...
    51.455157,-2.585359;...
    51.455157,-2.585359; 
    ];    


txLocations = TxLocations(Txroute,noLoc); % TXs returns noLoc txsites
                                         %  locations specified by the 
                                         %  rows of R (R for route).
                                         
                                        
txs= txsite("Name","TX Site", ...
    "Latitude",txLocations(:,1), ...
    "Longitude",txLocations(:,2), ...  
    "AntennaHeight",1, ...
    "TransmitterPower",10^(txpower/10)/1000, ... % 1mW or 1d\phbm default is 5 Watts
    "TransmitterFrequency", f);

rxs = rxsite("Name",rxNames, ...
    "Latitude",rxLocations(:,1), ...
    "Longitude",rxLocations(:,2), ...
    "ReceiverSensitivity",-120,...
    "AntennaHeight",[95,105]);

show(rxs); show(txs);

%% Propagation model
rtpm = propagationModel("raytracing", ...
    "Method","sbr", ... %or image for exact locations but option limited for up to reflections.
    "MaxNumReflections",maxNoRef, ...
    "BuildingsMaterial","perfect-reflector", ...
    "TerrainMaterial","perfect-reflector");
rays = raytrace(txs,rxs,rtpm);

msg = "no rays to plot"
for t = 1:length(txs)
    show(txs(t))
    for r = 1:length(rxs)
        try
            plot(rays{t,r})
        catch 
            display(msg)
        end
    end
end

% % Signal strength at the receivers built-in function
% ss = []; % i^th row of SS stores the receives signal strength at the ith recevier 
% for t = 1:length(txs)
%     for r = 1:length(rxs)
%         ss(t,r) = sigstrength(rxs(r),txs(t),rtpm);
%     end
% end
% ssLinear = 10.^(ss/10);  % mWatts
% 
%% Signal strength at the receivers my function
ss = []; % i^th row of SS stores the receives signal strength at the ith recevier 
cumulPhaseRx = []; cumulAttenuation = []
for t = 1:length(txs)
    for r = 1:length(rxs)
        ray = rays{t,r};
        NoRef{r,t} = [ray.NumReflections];
        PL = [ray.PathLoss]; % pathloss specified for each ray
        PH = [ray.PhaseShift];
        summReal = 0; summIm = 0;
        if length(ray) ~= 0
           % rtchan{r} = comm.RayTracingChannel(rays{r},txs(t),rxs(r));
           % showProfile(rtchan{i})
            for k=1:length(PH)
                summReal = summReal + 10^(-PL(k)/10) * cos(PH(k));
                summIm   = summIm   + 10^(-PL(k)/10) * sin(PH(k));
            end
            complexGain(t,r)= complex(summReal, summIm);
            cumulPhaseRx(t,r) = angle(complexGain(t,r));
            cumulAtten = abs(complexGain(t,r));
            cumulAttenuation(t,r) = -10*log10(cumulAtten);
        end
            ss = txpower - cumulAttenuation;
            ssLinear = 10.^(ss/10);  % mWatts
    end
end

% %% Channel statistics
% ssMean = mean(ssLinear);
% ssVar = var(ssLinear);
% meanProd = mean(ssLinear(:,1).*ssLinear(:,2)); 
%% Channel correlation -  REVISED
R = corr(complexGain(:,1),complexGain(:,2))
Rsqrd = abs(R)^2
Renvelope = corr(abs(complexGain(:,1)),abs(complexGain(:,2)))
Rpower = corr(abs(complexGain(:,1)).^2,abs(complexGain(:,2)).^2)
dist = distance(rxs(1),rxs(2))/lambda

%% Write data to excel file
if write == "yes"
    newData = {fc/1e6, maxNumReflections, length(txs),...
                distrx, real(R),imag(R),...
               Rsqrd, Renvelope, Rpower,...
               modelType," comment go here "};
    s = xlsappend('createfile.xlsx',newData);
    assert(s == 1, "addition of new data on the spreadsheet failed")
    % open('office.xlsx')
end

        
        
            
            