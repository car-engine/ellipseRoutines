%% distBetweenLL
% Perry Hong
% 22 March 2023
%
% Evaluates the distance (m) between two LL coordinate points 
% INPUT: coordinates = [long lat]
% Assumes WGS84 spheroid
% t = latitude in degrees
% 1 deg latitude = 111132.92 - 559.82 cos(2t) + 1.175 cos(4t) - 0.0023 cos(6t)
% 1 deg longitude = 111412.84 cos (t) - 93.5 cos (3t) + 0.118 cos (5t)

%% Begin function
function distOut = distBetweenLL(coordinates1, coordinates2)

    long1 = coordinates1(1);
    lat1 = coordinates1(2);
    long2 = coordinates2(1);
    lat2 = coordinates2(2);
    
    averageLatRad = deg2rad((lat1 + lat2)/2);

    latScale = 111132.92 - 559.82*cos(2*averageLatRad) + 1.175*cos(4*averageLatRad) - 0.0023*cos(6*averageLatRad);
    longScale = 111412.84*cos(averageLatRad) - 93.5*cos(3*averageLatRad) + 0.118*cos(5*averageLatRad);
    
    latDist = abs(lat2 - lat1)*latScale;
    longDist = abs(long2 - long1)*longScale;
    
    distOut = sqrt(latDist^2 + longDist^2);
    
end
