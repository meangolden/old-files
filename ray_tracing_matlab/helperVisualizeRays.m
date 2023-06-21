function helperVisualizeRays(rays)
%HELPERVISUALIZERAYS Visualize rays in a 3-D map
%   helperVisualizeRays(RAYS) visualizes the rays represented by the
%   comm.Ray objects, RAYS, in a 3-D map. The rays are colored according to
%   their path loss values.

%   Copyright 2020 The MathWorks, Inc. 

% Grid the path loss range with 5 dB separation
PL = [rays.PathLoss];
minPL = floor(min(PL)/10)*10;
maxPL = ceil(max(PL)/10)*10;
numColors = (maxPL - minPL)/5;
cm = colormap(jet(numColors+1));

% Plot rays
for i = 1:length(rays)
    if rays(i).LineOfSight    
        propPath = [rays(i).TransmitterLocation, rays(i).ReceiverLocation];
    else
        propPath = [rays(i).TransmitterLocation, ...
            rays(i).ReflectionLocations, rays(i).ReceiverLocation];
    end
    
    cmIdx = find((PL(i) - (maxPL:-5:minPL)) > 0, 1);
    line(propPath(1,:), propPath(2,:), propPath(3,:), ...
        'Color', cm(cmIdx,:), 'LineWidth', .8);
end

% Add color bar
cb = colorbar; 
cb.Label.String = 'Path Loss (dB)';
cbLim = cb.Limits;
cb.Ticks = cbLim(1) + diff(cbLim)/(2*numColors) + ...
    (0:numColors-1)*diff(cbLim)/numColors;
cb.TickLabels = maxPL-5:-5:minPL;

% [EOF]