function [ wb, t ] = gyroscope( w ) 
%gyroscope      Gyroscope angular velocity reading
%
%   [wb] = gyroscope(w) will produce the body-coordinate angular velocity
%   vector 'wb' by adding gyro bias and gaussian noise to the model angular
%   velocity vector 'w'. 
%
%   [wb, t] = gyroscope(w) will also output the orientation change due to
%   the angular velocity vector 'wb' in 't'. The output 't' will be
%   calculated from 'wb' to account for gyroscopic drift in the sensor
%   reading.
%
%   Class support for input w
%     float: double, single

if nargin < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

% Define sensor bias and noise parameters
bias = []; accuracy = []; sigma = [];
wb = (accuracy/sigma)*randn(3,1)+bias+w;

end