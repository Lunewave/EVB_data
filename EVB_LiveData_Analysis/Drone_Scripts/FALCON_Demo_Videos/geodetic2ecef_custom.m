function [x, y, z] = geodetic2ecef_custom(lat, lon, h)
%GEODETIC2ECEF_CUSTOM Convert geodetic coordinates to ECEF (Earth-Centered, Earth-Fixed)
%
%   [x, y, z] = geodetic2ecef_custom(lat, lon, h)
%
%   Inputs:
%     lat - Latitude in degrees
%     lon - Longitude in degrees
%     h   - Height above ellipsoid in meters
%
%   Outputs:
%     x, y, z - ECEF coordinates in meters

    % WGS84 ellipsoid constants
    a = 6378137.0;                % semi-major axis [m]
    f = 1 / 298.257223563;        % flattening
    e2 = f * (2 - f);             % eccentricity squared

    % Convert angles to radians
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);

    % Prime vertical radius of curvature
    N = a ./ sqrt(1 - e2 * sin(lat_rad).^2);

    % ECEF coordinates
    x = (N + h) .* cos(lat_rad) .* cos(lon_rad);
    y = (N + h) .* cos(lat_rad) .* sin(lon_rad);
    z = (N .* (1 - e2) + h) .* sin(lat_rad);
end
