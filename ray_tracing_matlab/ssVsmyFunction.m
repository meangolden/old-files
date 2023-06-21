%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: MATLAB Propagation Simulation
% Author: [Your Name]
% Date: [Date]

% Description:
% This MATLAB code simulates propagation of radio waves in a specified
% environment using ray tracing. It calculates the signal strength received
% at a receiver from multiple transmitter locations and compares it with
% the built-in signal strength calculations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Parameters
clear all; close all;

% Initialize variables
TRANSMITTER_HEIGHT = 10;      % Transmitter antenna height (m)
RECEIVER_HEIGHT = 1.5;        % Receiver antenna height (m)
RECEIVER_SENSITIVITY = -80;   % Receiver sensitivity (dBm)

maxNoRef = 4;
txpower = 0;   % in dBm
f = 400e6;

%% Map and propagation model
mapFileName = "bristol_campus.osm";
viewer = siteviewer("Buildings", mapFileName, "Basemap", "topographic");

rtpm = propagationModel("raytracing", ...
    "Method", "sbr", ...
    "MaxNumReflections", maxNoRef, ...
    "BuildingsMaterial", "perfect-reflector", ...
    "TerrainMaterial", "perfect-reflector");

%% Transmitter's locations
txLocations = [...
    51.456815 -2.603504; ...
    51.456828 -2.603517; ...
    51.456841 -2.603520; ...
    51.456855 -2.603531; ...
    51.456874 -2.603545; ...
    51.456890 -2.603554; ...
     51.456900 -2.603565; ...
    51.456926 -2.603581; ...
    51.456944 -2.603593];

%% Receiver's properties
rx = rxsite("Name", "Bob", ...
    "Latitude", 51.457377, ...
    "Longitude", -2.602006, ...
    "ReceiverSensitivity", RECEIVER_SENSITIVITY, ...
    "AntennaHeight", RECEIVER_HEIGHT);

%% Main
BuiltinSS = zeros(size(txLocations, 1), 1);
MySS = zeros(size(txLocations, 1), 1);

for t = 1:size(txLocations, 1)
    txLatitude = txLocations(t, 1);
    txLongitude = txLocations(t, 2);
    
    tx = txsite("Name", "Alice", ...
        "Latitude", txLatitude, ...
        "Longitude", txLongitude, ...
        "AntennaHeight", TRANSMITTER_HEIGHT, ...
        "TransmitterPower", 10^(txpower/10)/1000, ...
        "TransmitterFrequency", f);
    
    ray = raytrace(tx, rx, rtpm);
    PL = pathloss(rtpm, rx, tx);
    
    PH = ray.PhaseShift;
    cumulAttenuation = -10 * log10(abs(sum(10.^(-PL/10) .* exp(1i * PH))));
    MySS(t) = txpower - cumulAttenuation;
    BuiltinSS(t) = sigstrength(rx, tx, rtpm);
end

%% Display on command window
cumulAttenuation
pl = PL
ss = sigstrength(rx, tx, rtpm);
disp("Rtxpower - cumulAttenuation: " + (txpower - cumulAttenuation) + " dB

