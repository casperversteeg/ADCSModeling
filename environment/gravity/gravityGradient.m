function Tgg = gravityGradient(PosN, A, J)
%gravityGradient    Returns gravity gradient torque on satellite
%
%   Tgg = gravityGradient(PosN, A, J) returns the gravity gradient
%   disturbance torque on a satellite. The vector PosN contains the
%   inertial coordinates of the satellite position, the matrix A is the
%   local attitude matrix of the satellite relative to the inertial frame,
%   and the matrix J contains the satellite moments of inertia. 
%
%   The satellite inertia can be reported in either principal diagonal, or
%   as off-axis inertia. This will be corrected with an eigenvalue function
%   and sorted for computation.
%
%   Class support for inputs PosN, A, J
%     float: double, single

if nargin < 3
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

% Express attitude in terms of roll-pitch-yaw by 1-2-3 Euler angles
RPY = attitudeMatrix2RPY123(A);
theta = RPY(2); psi = RPY(3);
% Compute principal inertias J1, J2 and J3
[J, ~] = eig(J);
J = diag(J);
J = sort(J, 'desc');
% Gravity constant
mu = 3.986004418e14;
% Gravity Gradient torque
Tgg = 3*mu/norm(PosN)*[(J(3)-J(2))*cosd(theta)^2*cosd(psi)*sind(psi);
    (J(3)-J(1))*cosd(theta)*sind(theta)*cosd(psi);
    (J(1)-J(2))*cosd(theta)*sind(theta)*sind(psi)];
Tgg = (0.9 + (1.1-0.9)*randn(3,1)).*Tgg;
end