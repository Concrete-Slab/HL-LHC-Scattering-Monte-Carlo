function P2 = lorentz(P1,betaFrame)
%LORENTZ Lorentz transform of a four-momentum along the z axis
%   Detailed explanation goes here
arguments(Input)
    P1(4,1) double          % Four-momentum in initial frame (GeV/c)
    betaFrame(1,1) double   % beta factor of the second frame rel to initial frame
end
arguments(Output)
    P2(4,1) double          % Four-momentum in second frame (GeV/c)
end
    gammaFrame = 1/sqrt(1-betaFrame^2);
    T = [gammaFrame 0 0 -gammaFrame*betaFrame; 0 1 0 0; 0 0 1 0; -gammaFrame*betaFrame 0 0 gammaFrame];
    P2 = T*P1;
end