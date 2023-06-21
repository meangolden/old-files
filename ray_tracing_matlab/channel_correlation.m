clear all; close all;
%% Main Parameters
maxNoRef= 8            ;
txpower= 0;  % in dbm
f = 1e9; % in Hz
noLoc = 20;

%% Map
mapFileName = "bristol_campus.osm";
viewer = siteviewer("Buildings",mapFileName ,"Basemap","topographic");

%% Tranceivers
rxNames = {'Bob', 'Eve'};
assert(length(rxNames)==2,"This code is programmed for two receivers only")
Txroute  = [...
            51.457674,-2.600435
            51.458335,-2.600397; ... % This the route followed by the tx;
            51.458693,-2.599845; ...
            51.459053,-2.600274; ...
            51.459500,-2.600927; ...
            51.459234,-2.601786; ...
            51.458940,-2.602697; ...
            51.458722,-2.603423; ...
            51.458567,-2.603784 ;... 
            51.458301,-2.603750 ;...
            51.458019,-2.603916;...
            51.457656,-2.603872;...
            51.457319,-2.603747;...
            51.456941,-2.603565;...
            ]; %  51.458019,-2.603916
        
rxLocations =[...
    51.457377 -2.602006 ;...
    51.457389 -2.602006  ];    %  51.457362 -2.601989

% rxLocations =[...
%     51.458272,-2.601511 ;...
%     51.458300,-2.601371  ];    %  51.457362 -2.601989

% rxLocations =[...
%     51.458272,-2.601511 ;...
%     51.456245,-2.602275  ];    %  51.457362 -2.601989

txLocations = TxLocations(Txroute,noLoc); % TXs returns noLoc txsites
                                         %  locations specified by the 
                                         %  rows of R (R for route).
                                         
                                        
txs= txsite("Name","TX Site", ...
    "Latitude",txLocations(:,1), ...
    "Longitude",txLocations(:,2), ...  
    "AntennaHeight",1.5, ...
    "TransmitterPower",10^(txpower/10)/1000, ... % 1mW or 1d\phbm default is 5 Watts
    "TransmitterFrequency", f);

rxs = rxsite("Name",rxNames, ...
    "Latitude",rxLocations(:,1), ...
    "Longitude",rxLocations(:,2), ...
    "ReceiverSensitivity",-90,...
    "AntennaHeight",1.5);

show(rxs); 

%% Propagation model
rtpm = propagationModel("raytracing", ...
    "Method","sbr", ...
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

%% Channel statistics
ssMean = mean(ssLinear);
ssVar = var(ssLinear);
meanProd = mean(ssLinear(:,1).*ssLinear(:,2)); % mean of the pairwise
                                               % channel product (E(h1*h2))
% %% Channel Correlation IT IS WRONG./ CHECK NEXT SECTION
% numer = meanProd - ssMean(1)*ssMean(2);
% % denom = sqrt(ssVar(1)*ssVar(2)) % It doesn't work well for a small sample
%                               %  perhaps it has to do with sth about
%                               %  unbiased vs biased.
% V = mean(ssLinear.^2 - mean(ssLinear).^2); % Tested. Seems to perform well!
% denom = sqrt(prod(V));
% chCorr = numer / denom

%% Channel correlation -  REVISED
R = corr(complexGain(:,1),complexGain(:,2))
dist = distance(rxs(1),rxs(2))



        
        
            
            