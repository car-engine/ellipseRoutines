%% checkInEllipse
% 15 March 2023
% Perry Hong
%
% Checks if (pointLong, pointLat) is within an error ellipse
% Ellipse equation: (((x-h)cos(t) + (y-k)sin(t))/a)^2 + (((x-h)sin(t) - (y-k)cos(t))/b)^2 = 1
% ellipseAngle (degrees) is measured from the x-axis, anticlockwise
% ellipseLong, ellipseLat, ellipseMajor, ellipseMinor, pointLong and pointLat are all in degrees 

%% Begin function

function inEllipse = checkInEllipse(ellipseLong, ellipseLat, ellipseMajor, ellipseMinor, ellipseAngle, pointLong, pointLat)
    
    arguments
    
        ellipseLong (:,1) double
        ellipseLat (:,1) double
        ellipseMajor (:,1) double
        ellipseMinor (:,1) double
        ellipseAngle (:,1) double
        pointLong (:,1) double
        pointLat (:,1) double

    end

    if size(pointLong,1) ~= 1 && (size(pointLong,1) ~= size(ellipseLong,1))
        error('pointLong and pointLat must either exactly one row, or the same number of rows as the ellipse inputs')
    end

    ellipseAngle = deg2rad(ellipseAngle);

    point = (((pointLong-ellipseLong).*cos(ellipseAngle) + (pointLat-ellipseLat).*sin(ellipseAngle))./ellipseMajor).^2 + (((pointLong-ellipseLong).*sin(ellipseAngle) - (pointLat-ellipseLat).*cos(ellipseAngle))./ellipseMinor).^2;

    inEllipse = (point <= 1);

end
