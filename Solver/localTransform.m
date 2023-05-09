function T = localTransform(c)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% c is the z direction of the local frame
c = c./norm(c);
globalDirs = eye(3);
if c'*globalDirs(:,1) ~= 1
    % c is not parallel to global x axis
    % first orthogonal vector is cross product of c with the global x axis
    b = cross(c,globalDirs(:,1));
else
    % c is parallel to global x axis
    % first orthogonal vector is cross of c with global y axis
    b = cross(c,globalDirs(:,2));
end
b = b./norm(b);
a = cross(b,c);
T = [a,b,c];
end