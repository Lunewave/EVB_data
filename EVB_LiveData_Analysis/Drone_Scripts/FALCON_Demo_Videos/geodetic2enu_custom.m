function [east, north, up] = geodetic2enu_custom(lat, lon, alt, ref_lat, ref_lon, ref_alt)
% Convert geodetic coordinates (lat, lon, alt) to local ENU coordinates
% relative to a reference point (ref_lat, ref_lon, ref_alt).
%
% Inputs:
%   lat, lon, alt   - vectors of latitude (deg), longitude (deg), altitude (m)
%   ref_lat, ref_lon, ref_alt - reference geodetic point (deg, deg, m)
%
% Outputs:
%   east, north, up - vectors of local ENU coordinates in meters

% Convert reference point to ECEF
[ref_x, ref_y, ref_z] = geodetic2ecef_custom(ref_lat, ref_lon, ref_alt);

% Convert input points to ECEF
[x, y, z] = geodetic2ecef_custom(lat, lon, alt);

% Compute delta vectors
dx = x - ref_x;
dy = y - ref_y;
dz = z - ref_z;

% Convert degrees to radians
phi = deg2rad(ref_lat);
lambda = deg2rad(ref_lon);

% Rotation matrix ECEF to ENU
R = [-sin(lambda),             cos(lambda),              0;
     -sin(phi)*cos(lambda),  -sin(phi)*sin(lambda),    cos(phi);
      cos(phi)*cos(lambda),   cos(phi)*sin(lambda),    sin(phi)];

% Preallocate outputs
nPoints = length(dx);
east = zeros(nPoints,1);
north = zeros(nPoints,1);
up = zeros(nPoints,1);

% Rotate each delta vector into ENU
for i = 1:nPoints
    enu = R * [dx(i); dy(i); dz(i)];
    east(i) = enu(1);
    north(i) = enu(2);
    up(i) = enu(3);
end
end
