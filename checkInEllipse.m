%% checkInEllipse
% Checks if (currentLong, currentLat) is within an error ellipse
% Ellipse equation: (((x-h)cos(t) + (y-k)sin(t))/a)^2 + (((x-h)sin(t) - (y-k)cos(t))/b)^2 = 1
% ellipseAngle (degrees) is measured from the x-axis, anticlockwise
% ellipseMajor and ellipseMinor are in degrees

%% Begin function

function inEllipse = checkInEllipse(ellipseLat, ellipseLong, ellipseMajor, ellipseMinor, ellipseAngle, currentLat, currentLong)
                   
    ellipseAngle = deg2rad(ellipseAngle);
    
    point = (((currentLong-ellipseLong)*cos(ellipseAngle) + (currentLat-ellipseLat)*sin(ellipseAngle))/ellipseMajor)^2 + (((currentLong-ellipseLong)*sin(ellipseAngle) - (currentLat-ellipseLat)*cos(ellipseAngle))/ellipseMinor)^2;

    inEllipse = (point <= 1);

end
