function [ v_cross ] = crossMatrix(v)
%crossMatrix    Find the cross matrix of a vector v
%
%   v_cross = crossMatrix(v) will return the cross matrix of a vector 'v',
%   where the cross matrix is defined as v_cross = [0 -v(3) v(2); v(3) 0
%   -v(1); -v(2) v(1) 0];
%
%   Class support for input v
%     float: double, single

if nargin < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

v_cross = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end