function [ q, b ] = A2q( A )
%A2q    Finds the quaternion representation of an attitude matrix
%
%   q = A2q(A) returns the rotation quaternion 'q' associated with a
%   rotation matrix 'A'. 
%
%   Class support for inputs A
%     float: single, double

if nargin < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

q = zeros(4,1);
q(4) = sqrt(1+trace(A))/2;
q(1) = 1/4/q(4)*(A(2,3)-A(3,2));
q(2) = 1/4/q(4)*(A(3,1)-A(1,3));
q(3) = 1/4/q(4)*(A(1,2)-A(2,1));

end