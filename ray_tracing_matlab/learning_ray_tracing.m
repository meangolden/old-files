mapFileName = "conferenceroom.stl";
fc = 915e6;
lambda = physconst('lightspeed')/fc;

% txAntenna = dipole; % takes ages!(never let it run till the end)
% rxAntenna = dipole;
% txAntenna = horn; % takes ages!(never let it run till the end)
% rxAntenna = horn;
txAntenna = arrayConfig("Size",[1 4],'ElementSpacing',2*lambda);
rxAntenna = arrayConfig("Size",[1 4],'ElementSpacing',lambda);

helperViewArray(txAntenna);
helperViewArray(rxAntenna);

tx = txsite("cartesian", ...
    "Antenna",txAntenna, ...
    "AntennaPosition",[-1.46; -1.42; 0.1], ...
    'TransmitterFrequency',915e6);

rx1 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[-.15; .3; .85], ...
    "AntennaAngle",[0;90]);

rx2 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[.18; .3; .85], ...
    "AntennaAngle",[0;90]);

rx3 = rxsite("cartesian", ...
    "Antenna",rxAntenna, ...
    "AntennaPosition",[.28; .3; .85], ...
    "AntennaAngle",[0;90]);

helperVisualizeScenario(mapFileName,tx,[rx1,rx2,rx3])
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...   # sbr:shooting and bouncing rays
    "AngularSeparation","low", ...
    "MaxNumReflections",2, ...
    "SurfaceMaterial","wood");

PL = cell(1); PH = []; NoRef = cell(1); PropDist = []; % initialisation
rays = raytrace(tx,[rx1,rx2,rx3],pm,'Map',mapFileName);

for i=[1,2,3]
    ray = rays{i};
    NoRef{i} = [ray.NumReflections];
    PropDist = [ray.PropagationDistance];
    PL{i} = [ray.PathLoss];
    PH{i} = [ray.PhaseShift]
    helperVisualizeRays(ray);
end


summReal = [0,0,0]; summIm = [0,0,0]; cumulPhaseRx=[]; cumulAttenuation=[];
for j = [1,2,3]
    for i=1:length(PH{j})
       summReal(j) = summReal(j) + 10^(-PL{j}(i)/10)*cos(PH{j}(i));
       summIm(j) = summIm(j) + 10^(-PL{j}(i)/10)*sin(PH{j}(i));
    end
    cumulPhaseRx(j) = mod(atan(summIm(j)/summReal(j)),2*pi);
    cumulAttenuation(j) = (summIm(j)^2 + summReal(j)^2)^0.5;
    cumulAttenuation(j) = -10*log10(cumulAttenuation(j));
end
cumulPhaseRx
cumulAttenuation

%continue from here. Following gives errors
rtchan = comm.RayTracingChannel(rays{1},tx,rx1)
showProfile(rtchan)