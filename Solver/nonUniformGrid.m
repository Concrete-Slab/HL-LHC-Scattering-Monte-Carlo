function nonUniform = nonUniformGrid(uniformgrid,densityFactor,func)
%NONUNIFORMGRID Summary of this function goes here
%   Detailed explanation goes here
% func must be monotonically increasing
xmax = max(uniformgrid);
xmin = min(uniformgrid);
transformer = @(x) (func(x.*densityFactor))./(func(densityFactor)); % ensures function evaluated at 1 is always 1
nonUniform = zeros(1,length(uniformgrid));
nonUniform(uniformgrid>0) = transformer(uniformgrid(uniformgrid>0)./xmax)*min(xmax,xmax-xmin)+max(xmin,0);
nonUniform(uniformgrid<0) = transformer(uniformgrid(uniformgrid<0)./xmin)*max(xmin,xmin-xmax)+min(xmax,0);
end

