function [ sF ] = fineSunSensor( s, A, B )
% fineSunSensor    Produces fine sun sensor out given a model sun vector
%   
%   [s] = findSunSensor(s, A, B) will return the output of a fine sun  
%   sensor, given an inertial sun vector 's', and a current attitude 
%   matrix 'A' to transform into a body-coordinate sun vector. The matrix 
%   B is the body-centered tranformation matrix for the position of the 
%   sun sensor.
%
%   For the matrix B, the third column should contain the coordinates for
%   the vector perpendicular to the sensor of the FSS, defined in the body
%   coordinate frame. The coordinates defined by the matrix B should be
%   orthogonal.
%
%   Class support for inputs s, A, B
%     float: double, single

if nargin < 3
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

% Transform sun vector from inertial to body to fine sun sensor
sF = B'*A*s;
% Check if FSS "sees" sun (FOV +/- 60 deg):
if acosd(dot(sF, [0; 0; 1])) > 180
    sF = zeros(3,1);
else
    % Transform sun vector into sun sensor coordinate system
    % sB = B*sB;
    % Compute sun coordinate angles:
    alpha = acosd(dot(sF, [1; 0; 0]));
    beta = acosd(dot(sF, [0; 1; 0]));
    % Add normally distributed Gaussian noise to angle readings (zero mean,
    % 0.5 degree 3-sigma CI)
    alpha = (0.5/3)*randn(1) + alpha;
    beta = (0.5/3)*randn(1) + beta;
    % Find direction cosines. Gamma dependent on error in alpha and beta.
    calpha = cosd(alpha); cbeta = cosd(beta);
    cgamma = sqrt(1-calpha^2-cbeta^2);
    % Compute sun vector measurement in body-coordinates
    FSS = [calpha; cbeta; cgamma];
    sF = B*FSS;
end
end