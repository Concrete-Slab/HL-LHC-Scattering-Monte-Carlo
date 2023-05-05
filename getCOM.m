function [betaCOM,minvariant] = getCOM(p1,E1,p2,E2)
%getCOM Computes the center of mass beta factor for a 3D system of two
%particles
arguments(Input)
    p1(1,1) double % particle 1 momentum (GeV/c)
    E1(1,1) double % particle 1 energy (GeV)
    p2(1,1) double % particle momentum (GeV/c)
    E2(1,1) double % particle 2 energy (GeV)
end
arguments(Output)
    betaCOM(1,1) double     % centre of mass beta factor
    minvariant(1,1) double  % invariant mass of system
end
pTot = sym(p1+p2);
ETot = sym(E1+E2);
mTot = sqrt(ETot^2-pTot^2);

betaCOM = double(sqrt((pTot./mTot)^2./(1+(pTot./mTot)^2)));
minvariant = double(mTot);
end