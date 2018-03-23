function [ T, F ] = aerodynDrag(VelN, A)
% aerodynDrag     Determines solar radiation pressure at given PosNition
%
%   [T, F] = aerodynDrag(PosN, VelN, A, epoch, PosOrder) returns the
%   torque vector 'T', and force vector 'F', in inertial coordinates, given
%   a velocity in an inertial reference frame (in Cartesian coordinates) 
%   specified by 'VelN' (1x3), and an attitude matrix 'A' (3x3).
%
%   Because of a lack of accurate aerodynamic disturbance models, the
%   magnitude of the force will be determined from the drag equation, and
%   then multiplied by a random function that causes variance in the
%   disturbance torque.
%
%   By default, the function solarRadPress will convert coordinate system
%   back to Cartesian, so the output disturbance torque will be in
%   inertial, Cartesian coordinates.
%
%   Class support for inputs VelN, A
%     float: double, single

% Check input arguments. Set defaults if variables not defined.
if nargin < 2
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

% Make sure VelN is a column vector
if ~iscolumn(VelN)
    VelN = VelN';
end

% Atmospheric density at 400km altitude per the MSISE-90 Model of Earth's
% Upper Atmosphere (kg m^-3)
LowActivityDensity = 5.68e-13; HighActivityDensity = 5.04e-11;
% Assign local density a psuedorandom value on the interval of low solar
% activity density to hihgh solar activity density (kg m^-3). Can update
% this to resemble a Gaussian distribution.
rho = LowActivityDensity + (HighActivityDensity - LowActivityDensity)...
    *randn(1);
% Define momentum transfer coefficients
s_n = 0.5; s_t = 0.5;
% Surface areas of sides of CubeSat; Area vector
A_3U = 0.3*0.1; A_1U = 0.1*0.1; 
surf = [A_1U, A_3U, A_3U]; surf = [surf, surf];
% Create velocity unit vector and velocity magnitude in body coordinates
V = norm(VelN); e_v = A*(VelN/V);
% Create check for frontal areas of the satellite (c > 0)
e_vBdot = [e_v, e_v, e_v, e_v, e_v, e_v]';
Bni = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
c = dot(Bni, e_vBdot, 2);
% Calculate aerodynamic drag force
F = 0;
for i = 1:6
    if c(i) > 0
        F = F + (-rho*V^2*((2-s_n-s_t)*(dot(e_v, Bni(i,:)))^2*Bni(i,:)'+...
            s_t*(dot(e_v, Bni(i,:)))*e_v)*surf(i));
    end
end
% Add 10% pseudorandom variation to the force vector due to uncertainty
% in solar pressure value, optical behavior, and moment arm to center
% of pressure
F = (0.9 + (1.1-0.9)*randn(3,1)).*F;
% Define center of pressure moment arm r, and the aerodynamic disturbance 
% torque T as the cross product rxF, all in body coordinates
r = [0.01; 0.01; 0.01];
T = cross(r,F);
end