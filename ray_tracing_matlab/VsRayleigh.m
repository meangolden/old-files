function [outputArg1,outputArg2] = VsRayleigh(ssLinear)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sz = size(ssLinear)
msg1 = "Only one channel at the time can be compared to a Rayleigh channel"
msg2 = "Signal strength must be converted to mW (linear)"
assert(size(1,2)== 1, msg1)
assert(ssLinear > -1);
ssMean = mean(ssLinear);
ssVar = var(ssLinear);


outputArg1 = gfd
outputArg2 = inputArg2;
end

