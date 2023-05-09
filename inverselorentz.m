function P1 = inverselorentz(P2,betaFrame)
%LORENTZ Lorentz transform of a four-momentum along the z axis
%   Detailed explanation goes here
arguments(Input)
    P2(4,1) double          % Four-momentum in second frame (GeV/c)
    betaFrame(1,1) double   % beta factor of the second frame rel to initial frame
end
arguments(Output)
    P1(4,1) double          % Four-momentum in initial frame (GeV/c)
end
    gammaFrame = 1/sqrt(1-betaFrame^2);
    Tinv = [gammaFrame 0 0 gammaFrame*betaFrame; 0 1 0 0; 0 0 1 0; gammaFrame*betaFrame 0 0 gammaFrame];
    P1 = Tinv*P2;
end