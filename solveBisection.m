function [g2,fg2] = solveBisection(f,g0,g1,abstol,maxIterations)
%SOLVEBISECTION Summary of this function goes here
%   Detailed explanation goes here
fg0 = f(g0);
fg1 = f(g1);
if fg0*fg1>0
    error("Guesses have same sign")
end
propagate = true;
i=0;
while propagate
    i=i+1;
    g2 = (g0+g1)/2;
    if abs(g1-g0)<abstol || i>=maxIterations
        propagate=false;
    else
        fg2 = f(g2);
        if fg0*fg2>0
            g0=g2;
            % change to fg0 = fg2!!!!!
            fg0 = fg2;
        else
            g1=g2;
        end
    end
    fprintf("Iteration %d, current guess is %.1f microns\n",i,g2*1e6)
end
end

