function I = getHaloIntensity(sigma,options)
%GETHALOINTENSITY Constant halo intensity, normalised or non-normalised
%   Given a gaussian core of rms width sigma, the halo intensity is 100,000
%   times less than the peak intensity in the core. If the halo is taken to
%   also contain 5% of the total energy, then the normalisation option will
%   account for this.
%   Returns the halo intensity in protons per metre per second
arguments(Input)
    sigma(1,1) double {mustBePositive} % beam core rms width, metres
    options.normalise(1,1) logical = 0;
end
arguments(Output)
    I(1,1) double % constant halo intensity
end
g0 = 1e-5 * 1/(2*pi*sigma^2);
if options.normalise
    Rcore = sigma*sqrt(-2*log(sigma^2*g0*2*pi));
    coreProbability = normcdf(Rcore,0,sigma);
    k = 0.95/coreProbability;
else
    k = 1;
end
I =  Consts.data.Nb/Consts.data.bunchSpacing * g0 * k;
end

