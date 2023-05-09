function [g2,fg2] = solveBisection(f,g0,g1,abstol,maxIterations)
%SOLVEBISECTION Summary of this function goes here
%   Detailed explanation goes here
tStart = tic;
fg0 = f(g0);
tElapsed = toc(tStart);
seconds = rem(tElapsed,60);
minutes = (tElapsed-seconds)./60;
fprintf(" - First guess is %.1f microns. Simulation time was %d mins %.5f secs\n",g0*1e6,minutes,seconds)
tStart = tic;
fg1 = f(g1);
tElapsed = toc(tStart);
seconds = rem(tElapsed,60);
minutes = (tElapsed-seconds)./60;
fprintf(" - Second guess is %.1f microns. Simulation time was %d mins %.5f secs\n",g1*1e6,minutes,seconds)
if fg0*fg1>0
    error("Guesses have same sign")
end
propagate = true;
i=0;
while propagate
    i=i+1;
%     g2 = 10^((log10(g0)+log10(g1))/2);
    % linear interpolation between upper and lower log values
    g2 = abs(g1-g0)/2;
    if abs(g1-g0)<abstol || i>=maxIterations
        propagate=false;
    else
        tStart = tic;
        fg2 = f(g2);
        tElapsed = toc(tStart);
        seconds = rem(tElapsed,60);
        minutes = (tElapsed-seconds)./60;
        fprintf(" - Guess %d is %.1f microns.  Simulation time was %d mins %.5f secs\n",i+2,g2*1e6,minutes,seconds)
        if fg0*fg2>0
            g0=g2;
            % change to fg0 = fg2!!!!!
            fg0 = fg2;
        else
            g1=g2;
        end
    end
end
end

