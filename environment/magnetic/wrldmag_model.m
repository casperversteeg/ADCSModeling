% Calculate Earth's magnetic field given:
%
% @param height     Scalar value (in meters) above Earth's surface
% @param lat        Scalar geodetic latitude (in degrees) where north
%                   latitude is positive and south latitude is negative
% @param lon        Scalar geodetic longtitude (in degrees) where 
%                   longitude is positive and west longitude is negative
% @param dyear       Scalar decimal year (fractional to include parts of the
%                   year)
%
% @return           xyz: Magnetic field vector [Bx,By,Bz] in nanotesla
%                   h:   Horizontal intensity in nanotesla
%                   dec: Declination in degrees
%                   dip: Inclination in degrees
%                   f:   Total intensity in nanotesla

prompt = 'height';
height = input(prompt);
prompt = 'lat';
lat = input(prompt);
prompt = 'lon';
lon = input(prompt);
prompt = 'dyear';
dyear = input(prompt);

[xyz, h, dec, dip, f] = wrldmagm(height, lat, lon, dyear);


