%% Main parameters
maxNumReflections = 3;
mapFileName = "conferenceroom.stl";
fc = 100e6;
txPower= 0.001;  % in watts
% Coordinates fot Txsites
X = [-1.5:1:1.5];
Y = [-1.5:1:1.5];
Z = [2.5];
% txAntenna = dipole; % takes ages!(never let it run till the end)
% rxAntenna = dipole;
% txAntenna = horn; % takes ages!(never let it run till the end)
% rxAntenna = horn;
lambda = physconst('lightspeed')/fc;
txAntenna = arrayConfig("Size",[1 1],'ElementSpacing',lambda);
rxAntenna = arrayConfig("Size",[1 1],'ElementSpacing',lambda);
modelType = "non-LOS"; % Checked visualeScenario

%% Tranceivers

txs = TxsIndoors(X,Y,Z,fc,txAntenna,txPower);

rx1 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[0.2; .3; 0.55], ...
    "AntennaAngle",[0;90]);

rx2 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[.25; .3; 0.55], ...
    "AntennaAngle",[0;90]);

rxs = [rx1, rx2];
helperVisualizeScenario(mapFileName,txs,rxs)

% helperViewArray(txAntenna);
% helperViewArray(rxAntenna);

%% Propagation model
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...   # sbr:shooting and bouncing rays
    "AngularSeparation","low", ...
    "MaxNumReflections",maxNumReflections, ...
    "SurfaceMaterial","wood");

rays = raytrace(txs,rxs,pm,'Map',mapFileName);

%% Signal strength at the receivers my function + plot rays
T = length(txs);
R = length(rxs);
ssLinear = zeros(T,R); % i^th row of SS stores the receives signal strength at the ith recevier 
cumulPhaseRx = zeros(T,R); complexGain = zeros(T,R); 
phaseShift = zeros(T,R);  %It seems to be independent to the number of subpaths :( I don't think I can use it. Finding another way. 
myPhaseshift = zeros(T,R);
msg = "no rays to plot";

for t = 1:T
    for r = 1:R
        ray = rays{t,r};
        PL = [ray.PathLoss]; % pathloss specified for each ray
        PH = [ray.PhaseShift];
        summReal = 0; summIm = 0;
        if ~isempty(ray)
           % helperVisualizeRays(rays{t,r}) % plot ray
           % rtchan{r} = comm.RayTracingChannel(rays{r},txs(t),rxs(r));
           % showProfile(rtchan{i})     
            
%            
            for k=1:length(PH)
                summReal = summReal + 10^(-PL(k)/10) * cos(PH(k));
                summIm   = summIm   + 10^(-PL(k)/10) * sin(PH(k));
                myPhaseshift(t,r) = myPhaseshift(t,r) + ...
                    ray(k).PropagationDistance;
                phaseShift(t,r) = phaseShift(t,r) + ray(k).PhaseShift;
            end
        end
            complexGain(t,r)= complex(summReal, summIm);
            cumulPhaseRx(t,r) = angle(complexGain(t,r));
            ssLinear(t,r) = abs(complexGain(t,r));
            ssDB = -10*log10(ssLinear);
    end
end

%% How rich is the multipath environment?
s = 0;
for k = 1:length(txs)
    s = s + numel(rays{k});
end
averageNoPaths = s/length(txs); % The average number of paths at receiver 
                                 % (based on the rays reaching rx1).
%% Results
myPhaseShift = (2*pi/lambda).*myPhaseshift;
R = corr(complexGain(:,1),complexGain(:,2))
Rshift = mean(cos(phaseShift(:,1) - phaseShift(:,2)))
Rmyshift = mean(cos(myPhaseShift(:,1) - myPhaseShift(:,2))) %gives the same result as Rshift. Cool
dist = distance(rxs(1),rxs(2))/lambda % normalised to wavelength
% %% Channel statistics
% if length(txs) == 1
%     ssMean = ssLinear;
% else
%     ssMean = mean(ssLinear);
% end
% meanProd = mean(ssLinear(:,1).*ssLinear(:,2)); % mean of the pairwise
%                                           % channel product (E(h1*h2))
% V = mean(ssLinear.^2 - mean(ssLinear).^2); % Tested. Seems to perform well!
% %% Channel Correlation
% assert(length(ssMean)==2, "exactly two rxs are required")
% numer = meanProd - ssMean(1)*ssMean(2);
% denom = sqrt(prod(V));
% chCorr = numer / denom

%% Write data to excel file

newData = {10*log10(txPower*1e3), fc/1e6, maxNumReflections, ...
           averageNoPaths,length(txs), dist, R, Rshift, modelType};
xlsappend('indoors.xlsx',newData);
%open('indoors.xlsx')


%% OK What I learnt today (23/07/2021)
% the build in fuction of corr(abs(complex gain at receiver 1)*abs(complex
% gain at receiver 2)) is the same as my function, which proves the my
% function works well. Still, things do not look well... though..... I am
% sorry, but progress hasn't beem made as you thought....What I need to do next:
% DO NO RELY on the phase information that the built in function gives you.
% It is random, and this could cause problems? I don't know.. What you
% have to do is to find the propagation distance of each ray. SO, find the
% phase by yourself. The location is not accurate, I heard you saying. Just
% experiment. And see if you can put the exact location methof in the
% indoor environemnt. The limitation is up to 2 reflections, I know...,
% Just investigate a bit more. One more day. Maximum 2 days. Good luck
% Chrys :)