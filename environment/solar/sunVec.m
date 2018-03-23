function [ s, Rs ] = sunVec(PosN, t, posOrder, epoch)
% sunVec    Produces a sun-pointing vector in satellite's body frame
%   
%   [s, Rs] = sunVec(PosN, t, posOrder, epoch) returns the sun vector s in
%   inertial coordinates given a position in space (Cartesian, sphiercal or
%   cylindrical coordinates) specified by 'PosN' (1x3), and a time 't'
%   since the epoch of the model. The output Rs is a scalar that gives the
%   instantaneous distance from the Earth to the Sun in AU.
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
%   By default, this model for the sun vector will use a spherical Earth
%   approximation to determine if the satellite is in the Earth's eclipse.
%   If the eclipse condition is true, the function will return a zero
%   vector.
%
%   Class support for inputs PosN, t, epoch
%     float: double, single
%   Class support for input posOrder
%     string: string

if nargin < 4
    epoch = [2015; 1; 1; 0; 0; 0];
    if nargin < 3
        posOrder = 'cart';
        if nargin < 2
            error(message('MATLAB:narginchk:notEnoughInputs'));
        end
    end
end

if strcmp(posOrder, 'sphere')
    PosN = sph2cart(PosN(1), PosN(2), PosN(3));
elseif strcmp(posOrder, 'cylin')
    PosN = pol2cart(PosN(1), PosN(2), PosN(3));
end

% Make epoch correction vector
if ~iscolumn(epoch)
    epoch = epoch';
end
if length(epoch) < 6
    for i = length(epoch)+1:6
        epoch(i) = 0;
    end
end

% Make time correction vector
SerialEpoch = datenum(epoch(1), epoch(2), epoch(3), epoch(4), epoch(5),...
    epoch(6));
tDays = t/24/3600;
ti = datevec(SerialEpoch + tDays);

% Correct current time for epoch and time in seconds since January 1, 2015
Y = ti(1); M = ti(2); D = ti(3); h = ti(4); m = ti(5); s = ti(6);
Tut1 = (JD(Y,M,D,h,m,s)-2451545)/36525;

% Source: Markley, F. Landis, Crassidis, John L. "Fundamentals of
% Spacecraft Attitude Determination and Control" Microcosm Press and
% Springer, 2014.
% Obliquity of Earth's rotational axis to orbital plane (degrees)
e = 23.439291-0.0130042*Tut1; 
% Radius of Earth (m)
R = 6.371e6;
% Mean longitude and mean anomaly of sun (degrees)
phi = 280.460 + 36000.771*Tut1;
M = 357.5277233 + 35999.05034*Tut1;
% Ecliptic longitude angle (degrees)
phi_e = phi + 1.91466671*sind(M)+0.019994643*sind(2*M);
% Sun-pointing vector in inertial coordinates
s = [cosd(phi_e); cosd(e)*sind(phi_e); sind(e)*sind(phi_e)];
% Compute instananeous distance from Earth to Sun
Rs = 1.000140612-0.016708617*cosd(M)-0.000139589*cosd(2*M);
% Find closest point to PosN on solar eclipse cylinder centerline
linePos = dot(PosN,s)*s;
% Check if the point distance is less than the radius of the Earth and lies
% in the eclipsed region. If it does, set sun vector to 0 (zero).
if (dist(linePos, PosN) < R) && (dot(s, PosN) < 0)
    s = zeros(3,1);
end

    function [dist] = dist(A, B)
        dist = sqrt((A(1)-B(1))^2+(A(2)-B(2))^2+(A(3)-B(3))^2);
    end
end