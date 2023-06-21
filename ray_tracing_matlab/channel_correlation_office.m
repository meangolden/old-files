
clear all;
for fc = 2e9%[0.1e9:0.1e9:8e9]
    clear AoA1 AoA2 
%% Main parameters
maxNumReflections = 5;  % do no go higher than two if progapation model is set to "image"                                                                                                                                                                                         ;
mapFileName = "office.stl";
fc
txPower= 0.001;  % in watts
material = "metal";
% "vacuum", "concrete",
%   "brick", "plasterboard", "wood", "glass", "ceiling-board", "chipboard",
%   "floorboard", "metal", "very-dry-ground", "medium-dry-ground" and "wet-
%   ground"
% Coordinates fot Txsites
X = [1:1:5]; % didn't start from 0 or they may be behind the bookshelves->zero channel coefficient
Y = [5:1:8];
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
    "AntennaPosition",[2; 1.5 ;.75], ...
    "AntennaAngle",[0;90]);

rx2 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[2.1; 1.5; .75], ...
    "AntennaAngle",[0;90]);

rxs = [rx1, rx2];
helperVisualizeScenario(mapFileName,txs,rxs)

% helperViewArray(txAntenna); helperViewArray(rxAntenna);

%% Propagation model
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...   # sbr:shooting and bouncing rays
    "AngularSeparation","low", ...
    "MaxNumReflections",maxNumReflections, ...
    "SurfaceMaterial",material); 
raysRx1 = raytrace(txs,rx1,pm,'Map',mapFileName);
raysRx2 = raytrace(txs,rx2,pm,'Map',mapFileName);

%% Complex channel gain
T = length(txs);
% initialisation
complexGainRx1 = zeros(T,1); complexGainRx2 = zeros(T,1); 
h1 = zeros(T,1); h2 = zeros(T,1);

% loops for every pair: transmitter-receiver
for t = 1:T
    for k = 1:numel(raysRx1{t})
%     pathLossRx1{t}{end + 1}  = 10^(- raysRx1{t}(k).PathLoss /10);
%     phaseShiftRx1{t}{end + 1}= raysRx1{t}(k).PhaseShift;
    complexGainRx1(t) = complexGainRx1(t) + ...
        10^(- raysRx1{t}(k).PathLoss/10) * exp(-raysRx1{t}(k).PhaseShift * j);
    end
    for k = 1:numel(raysRx2{t})
%     pathLossRx2{t}{end + 1}  = 10^(- raysRx2{t}(k).PathLoss / 10);
%     phaseShiftRx2{t}{end + 1}= raysRx2{t}(k).PhaseShift;
    complexGainRx2(t) = complexGainRx2(t) + ...
        10^(-raysRx2{t}(k).PathLoss /10) * exp(-raysRx2{t}(k).PhaseShift * j);
    end
    h1(t) = complexGainRx1(t);
    h2(t) = complexGainRx2(t);
end


%% more properties for a random selection of transmitter 
pathLossRx1 = []; pathLossRx2 = []; %in dB
phaseShiftRx1 = []; phaseShiftRx2 = [];
angleOfArrivalRx1 = []; angleOfArrivalRx2 = cell(T,1); %[azimuth, elevation]
noPaths = 0; % An indication of how rich the multipath is. 

randTx = randi([1,numel(txs)]); % for a random transmitter
while numel(raysRx1{randTx}) == 0
    randTx = randi([1,numel(txs)])
end
for k = 1:numel(raysRx1{randTx})
    pathLossRx1(k) = 10^(- raysRx1{randTx}(k).PathLoss /10); %path loss for each ray
    phaseShiftRx1(k)= raysRx1{randTx}(k).PhaseShift;
    AoA1(k,:) = (raysRx1{randTx}(k).AngleOfArrival);%*(pi/180); %comment out for rads
end
for k = 1:numel(raysRx2{randTx}) % for the same random transmitter
   pathLossRx2(k) = 10^(- raysRx2{randTx}(k).PathLoss / 10); %dB
   phaseShiftRx2(k)= raysRx2{randTx}(k).PhaseShift;
    noPaths = noPaths + numel(raysRx2{randTx}(k));
   AoA2(k,:) = (raysRx2{randTx}(k).AngleOfArrival);%*(pi/180); %comment out for rads
end
helperVisualizeScenario(mapFileName,txs(randTx),rxs)
ray1 = raysRx1{randTx}
helperVisualizeRays(ray1)
ray2 = raysRx2{randTx}
helperVisualizeRays(ray2)
%plot azimuth
figure
subplot(2,1,1)
polarscatter(AoA1(:,1),pathLossRx1)
title('Azimuth at Rx1')
subplot(2,1,2)
title('Azimuth at Rx2')
polarscatter(AoA2(:,1),pathLossRx2)
%plot elevation
figure
subplot(2,1,1)
polarscatter(AoA1(:,2),pathLossRx1)
subplot(2,1,2)
polarscatter(AoA2(:,2),pathLossRx2)

%% Channel correlation and Results
ComplexCorrelation = corr(h1,h2) 
ComplexCorrelationSqr = abs(ComplexCorrelation)^2
EnvelopeCorrelation = corr(abs(h1),abs(h2))
PowerCorrelation = corr(abs(h1).^2,abs(h2).^2)
distrx = distance(rxs(1),rxs(2))/lambda % normalised to wavelength
% averageAttenRx1 = 10*log10(mean(abs(h1))) %in dB
% averageAttenRx2 = 10*log10(mean(abs(h2))) %in dB

%% Write data to excel file
newData = {fc/1e6, maxNumReflections, noPaths,material, length(txs),...
            distrx, real(ComplexCorrelation),imag(ComplexCorrelation),...
           ComplexCorrelationSqr, EnvelopeCorrelation, PowerCorrelation,...
           modelType, "metal. rx2 next to rx1 x1:2, x2=2.1 "};
s = xlsappend('office.xlsx',newData);
assert(s == 1, "addition of new data on the spreadsheet failed")
% open('office.xlsx')

end


