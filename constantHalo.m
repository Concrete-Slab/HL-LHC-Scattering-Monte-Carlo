function [f,Ihalo,expectedProtons] = constantHalo(mu,sigma,cutoffFactor,xLimits,yLimits)
%HALODISTRIBUTION Summary of this function goes here
arguments                                   % all lengths in m
    mu(1,1) double                          % position of the centre of the beam
    sigma(1,1) double {mustBeNonnegative}   % rms width of beam core
    cutoffFactor(1,1) double {mustBePositive,mustBeLessThan(cutoffFactor,1)} % halo intensity relative to maximum intensity
    xLimits(2,:) double % x support region(s)
    yLimits(2,:) double % y support region(s)
end
% find total probability of gaussian up to cutoff value
g0 = 1/(2*pi*sigma^2);
ycut = g0*cutoffFactor;
rlim = sigma*sqrt(-2*log(sigma^2*ycut*2*pi));
gaussVolume = normcdf(rlim,0,sigma);
% find total probability elsewhere in the tube at constant pdf
constantVolume = ycut*pi*(Consts.TubeRadius^2-rlim^2);
if constantVolume<0
    error("Gaussian too wide to allow for constant halo within beam tube")
end
totalArea = constantVolume+gaussVolume;
constantVal = ycut./totalArea;
f = @(x) sqrt(constantVal);
% in a uniform distribution with independent supports, probability is the
% constant value * area of support
dy = sum(diff(yLimits));
dx = sum(diff(xLimits));
% totProb = constantVal*dx*dy
% expectedProtons = Ihalo*dx*dy*Consts.data.Nb./Consts.data.bunchSpacing;
Ihalo = g0*cutoffFactor*Consts.data.Nb./Consts.data.bunchSpacing;
% Ihalo = constantVal*Consts.data.Nb./Consts.data.bunchSpacing;


expectedProtons = Ihalo*dx*dy;
end

