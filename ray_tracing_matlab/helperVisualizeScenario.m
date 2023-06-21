function helperVisualizeScenario(mapFileName, txs, rxs)
%HELPERVISUALIZESCENARIO Visualize 3-D map and tx/rx sites
%   helperVisualizeScenario(MAPFILENAME, TXS, RXS) visualizes the 3-D map
%   in the STL file, MAPFILENAME, transmitter sites, TXS, in red, and
%   receiver site, RXS, in blue. 

%   Copyright 2020 The MathWorks, Inc. 

% Visualize the 3D map
if ~isa(mapFileName, 'triangulation')
    tri = stlread(mapFileName);
else
    tri = map;
end

% Visualize the 3D map from STL file
figure("Position", [360 360 600 600]);
trisurf(tri, ...
    'FaceAlpha', 0.3, ...
    'FaceColor', [.5, .5, .5], ...
    'EdgeColor', 'none');
view(60, 30);
hold on; axis equal; grid off;
xlabel('x'); ylabel('y'); zlabel('z');

% Plot edges
fe = featureEdges(tri,pi/20);
numEdges = size(fe, 1);
pts = tri.Points;
a = pts(fe(:,1),:); 
b = pts(fe(:,2),:); 
fePts = cat(1, reshape(a, 1, numEdges, 3), ...
    reshape(b, 1, numEdges, 3), nan(1, numEdges, 3));
fePts = reshape(fePts, [], 3);
plot3(fePts(:, 1), fePts(:, 2), fePts(:, 3), 'k', 'LineWidth', .5); 

% Visualize tx/rx sites
txPos = [txs.AntennaPosition];
rxPos = [rxs.AntennaPosition];
scatter3(txPos(1,:), txPos(2,:), txPos(3,:), 60, 'r', 'filled');
scatter3(rxPos(1,:), rxPos(2,:), rxPos(3,:), 60, 'b', 'filled');

end

% [EOF]