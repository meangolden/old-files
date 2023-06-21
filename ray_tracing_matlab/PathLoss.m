function [pl] = PathLoss(d,l)
% Returns the free path loss of a signal with wavelgnth lambda traveling d
% meters.
%   Variable d can be the propagation distance of a ray, in which case the
%   perfect reflection is considered. I.e., path loss due to difraction or
%   absorption are not taken into account. 
pl = (4*pi*d/l)^2
end

