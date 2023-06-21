function [txLocations] = TxLocations(R, noLoc)
%   TxLocations returns multiple txsites across the route specified by 
%   the rows of R.
%   TxLocations returns noLoc txsites at locations specified by the
%   rows of R. The transmit antenna is isotropic.
%   The geographical coordinates of R shall be given to be
%   roughly equidistant from oneanother. noBtw more geographical coordinates
%   are added in between each location at R. 
%   Thus, the output is noBtw*length(A)txsites.
%   Recall that the location accuracy in openstreetmap (and
%   therefore here) is 10cm which is provided by 6 dps at geographical
%   coordinates. As long as variable noBtw is a small power of ten,
%   it is expected that, the distance between any two txsites
%   is larger to that number.

t = 1:length(R');
pp = spline(t,R');

tInterp = linspace(1,length(R), noLoc);
txLocations = ppval(pp, tInterp);
txLocations = txLocations';
end
