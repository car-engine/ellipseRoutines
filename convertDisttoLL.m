%% convertDisttoLL
% Perry Hong
% 22 March 2023
%
% Converts a distance in SI units (m) to LL
% Assumes WGS84 spheroid
% t = latitude in degrees
% 1 deg latitude = 111132.92 - 559.82 cos(2t) + 1.175 cos(4t) - 0.0023 cos(6t)
% 1 deg longitude = 111412.84 cos (t) - 93.5 cos (3t) + 0.118 cos (5t)

%% Begin function
function [longLL, latLL] = convertDisttoLL(longDist, latDist, lat)

    latRad = deg2rad(lat);

    longScale = 111412.84*cos(latRad) - 93.5*cos(3*latRad) + 0.118*cos(5*latRad);
    latScale = 111132.92 - 559.82*cos(2*latRad) + 1.175*cos(4*latRad) - 0.0023*cos(6*latRad);

    longLL = longDist./longScale;
    latLL = latDist./latScale;
    

end
