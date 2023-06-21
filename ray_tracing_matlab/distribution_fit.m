
%% Main parameters
maxNumReflections = 2;  % do no go higher than two if progapation model is set to "image"                                                                                                                                                                                         ;
mapFileName = "office.stl";
txPower= 10;  % in watts
material = "metal";
fc = 2e9;
% Coordinates fot Txsites
X = (1:1:5); % didn't start from 0 or they may be behind the bookshelves->zero channel coefficient
Y = (5:1:8);
Z = (1:1:3);
lambda = physconst('lightspeed')/fc;
txAntenna = arrayConfig("Size",[1 1],'ElementSpacing',lambda);
rxAntenna = arrayConfig("Size",[1 1],'ElementSpacing',lambda);
modelType = "non-LOS"; % Checked visualeScenario

%% Tranceivers
txs = TxsIndoors(X,Y,Z,fc,txAntenna,txPower);
rx1 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[2; 1.5 ;.75], ...
    "AntennaAngle",[0;90]);

% helperVisualizeScenario(mapFileName,txs,rx1)
% helperViewArray(txAntenna); helperViewArray(rxAntenna);

%% Propagation model
clear AoA1
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...   # sbr:shooting and bouncing rays
    "AngularSeparation","low", ...
    "MaxNumReflections",maxNumReflections, ...
    "SurfaceMaterial",material); 
raysRx1 = raytrace(txs,rx1,pm,'Map',mapFileName);

%% Complex channel gain
T = length(txs);
% initialisation
complexGainRx1 = zeros(T,1);
h1 = zeros(T,1);

% loops for every pair: transmitter-receiver LOOP requires considerable time
for t = 1:T
    for k = 1:numel(raysRx1{t})
complexGainRx1(t) = complexGainRx1(t) + ...
        10^(- raysRx1{t}(k).PathLoss/10) * exp(-raysRx1{t}(k).PhaseShift * 1i);
    end
    h1(t) = complexGainRx1(t);
end

%% more properties for a random selection of transmitter 
pathLossRx1 = [];  %in dB
phaseShiftRx1 = []; 
angleOfArrivalRx1 = [];  %[azimuth, elevation]
noPaths = 0; % An indication of how rich the multipath is. 
randTx = randi([1,numel(txs)]); % for a random transmitter
while numel(raysRx1{randTx}) == 0
    randTx = randi([1,numel(txs)]); % in case there is no path for that txer
end

for k = 1:numel(raysRx1{randTx}) % for the same random transmitter
    pathLossRx1(k) = 10^(- raysRx1{randTx}(k).PathLoss / 10); %dB
   % phaseShiftRx2(k)= raysRx2{randTx}(k).PhaseShift;
    noPaths = noPaths + numel(raysRx1{randTx}(k))
    AoA1(k,:) = (raysRx1{randTx}(k).AngleOfArrival);%*(pi/180); %comment out for rads
end
helperVisualizeScenario(mapFileName,txs(randTx),rx1)
ray1 = raysRx1{randTx};
helperVisualizeRays(ray1)

%plot azimuth
figure
subplot(2,1,1)
polarscatter(AoA1(:,1),pathLossRx1)
title('Azimuth at Rx')
%plot elevation
subplot(2,1,2)
polarscatter(AoA1(:,2),pathLossRx1)
title('Elevation at Rx')

%% distribution fit



