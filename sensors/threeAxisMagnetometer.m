function [ B ] = threeAxisMagnetometer( Bn, A, TAM )
% threeAxisMagnetometer    Produces magnetometer sensor reading from model
%   
%   [B] = threeAxisMagnetometer(Bn, A, TAM) will return the output of a
%   magnetometer given a model magnetic field vector Bn in inertial 
%   coordinates. The output vector B will be in the body-coordinate system.
%   The matrix TAM (optional) is a rotation matrix that transforms the 
%   axes of the satellite to that of the magnetometer, if it is tilted 
%   off-axis. If TAM is not set, it will be assumed TAM = eye(3).
%
%   For the matrix TAM, the third column should contain the coordinates 
%   for the vector perpendicular to the sensor of the magnetometer, 
%   defined in the body coordinate frame. The coordinates defined by the 
%   matrix B should be orthogonal.
%
%   Class support for inputs s, A, B
%     float: double, single

if nargin < 3
    TAM = eye(3);
    if nargin < 2
        error(message('MATLAB:narginchk:notEnoughInputs'));
    end
end

% Transform magnetic field vector from inertial to body to magnetometer
Bbody = TAM'*A*Bn;
% Noise parameters: accuracy and associated confidence interval
accuracy = []; sigma = []; 
% Add noise to local magnetic field value for sensor output
B = (accuracy/sigma)*randn(3,1) + Bbody;

end