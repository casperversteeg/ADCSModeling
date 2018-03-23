function [T, F, s] = solarRadPress(PosN, A, t, posOrder, epoch)
% solarRadPress     Determines solar radiation pressure at given position
%
%   [T, F, s] = solarRadPress(PosN, A, t, posOrder, epoch) returns the
%   torque vector T, force vector F, and sun-pointing unit vector s, in 
%   inertial coordinates, given a position in space (Cartesian, spherical 
%   or cylindrical coordinates) specified by 'PosN' (1x3), and an attitude 
%   matrix 'A' (3x3) at a time 't' in seconds since the epoch. 
%
%   The input'epoch' determines the sun vector 's' in inertial coodinates, 
%   based on the location of Earth in its heliocentric orbit. The input
%   must be a vector in the form [Y; M; D; h; m; s]. If not set, the epoch
%   of the model will be 2015, January 1, 0:00:00. 
%   
%   The variable 'posOrder' is a setting for the 'PosN' input, to 
%   differentiate between the various coordinate types. It can take the
%   following values for the respective input types:
% 
%       'cart'   -   (default) Cartesian coordinate input for 'PosN'. 
%                    'PosN' should be a vector containing satellite 
%                    position in the form [x, y, z].
%       'sphere' -   Spherical coorindate input for 'PosN'. 'PosN' takes 
%                    the form [azimuth, elevation, altitude] to describe
%                    satellite position in orbit.
%       'cylin'  -   Cylindrical coordinate input for 'PosN'. For orbits
%                    described in cylindrical coordinates [azimuth, r, z].
%                    This is uncommon, but supported. 
%
%   By default, the function solarRadPress will convert coordinate system
%   back to Cartesian, so the output disturbance torque will be in
%   inertial, Cartesian coordinates.
%
%   Class support for inputs PosN, A, t, epoch
%     float: double, single
%   Class support for input posOrder
%     string: string

% Check input arguments. Set defaults if variables not defined.
if nargin < 5
    epoch = [2015; 1; 1; 0; 0; 0];
    if nargin < 4
        posOrder = 'cart';
        if nargin < 3
            error(message('MATLAB:narginchk:notEnoughInputs'));
        end
    end
end

if strcmp(posOrder, 'sphere')
    PosN = sph2cart(PosN(1), PosN(2), PosN(3));
elseif strcmp(posOrder, 'cylin')
    PosN = pol2cart(PosN(1), PosN(2), PosN(3));
end

% Specular reflection and diffusion coefficients for each of the six
% surfaces of a CubeSat. Surface property. (unitless in [0,1]) Current
% values are speculative.
Rspec = 0.6; Rspec = Rspec.*ones(1,6); 
Rdiff = 0.2; Rdiff = Rdiff.*ones(1,6);
% Surface areas of sides of CubeSat; Area vector
A_3U = 0.3*0.1; A_1U = 0.1*0.1; 
surf = [A_1U, A_3U, A_3U]; surf = [surf, surf];
% Get sun vector in inertial coordinates
[s, Rs] = sunVec(PosN, t, posOrder, epoch);
% Transform sun vector to body coordinates
sB = A*s;
% Create check for which sides of the satellite can be "seen" (c > 0)
sBdot = [sB, sB, sB, sB, sB, sB]';
Bni = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
c = dot(Bni, sBdot, 2);
% Solar radiation intensity variation, given a 5% pseudorandom uncertainty 
% fluctuation to account for solar activity and uncertainty in the
% parameter. Distance values must be in terms of AU.
Rs = s*Rs; R = Rs-(PosN./149597870700);
P = (1362/(299792458*norm(R)^2))*(0.95 + (1.05-0.95)*randn(1));
% Calculate solar pressure force
F = 0;
for i = 1:6
    if c(i) > 0
        F = F + (-P*surf(i)*c(i)*(2*(Rdiff(i)/3+Rspec(i)*c(i))*...
            (A'*Bni(i,:)')+(1-Rspec(i))*sB));
    end
end
% Add 10% pseudorandom variation to the force vector due to uncertainty
% in solar pressure value, optical behavior, and moment arm to center
% of pressure
F = (0.9 + (1.1-0.9)*randn(3,1)).*F;
% Define center of pressure moment arm r, and the radiation pressure 
% torque T as the cross product rxF, all in body coordinates
r = [0.01; 0.01; 0.01];
T = cross(r,F);
end