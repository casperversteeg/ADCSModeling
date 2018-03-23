function [ A ] = attitudeMatrix( b, n )
% attitudeMatrix    Computes a rotation matrix from two vectors
%
%   [A] = attitudeMatrix(b, n) returns the rotation matrix 'A' that
%   transforms the inertial-frame vector 'n' into the body-frame vector
%   'b'. 
%
%   Class support for inputs b, n
%     float: double, single

if nargin < 2
    error(message('MATLAB;narginchk:notEnoughInputs'));
end

if isrow(b)
    b = b';
end
if isrow(n)
    n = n';
end

v = cross(n,b);
s = norm(v);
c = dot(b,n);

A = eye(3) + crossMatrix(v) + (crossMatrix(v))^2*((1-c)/s^2);

end