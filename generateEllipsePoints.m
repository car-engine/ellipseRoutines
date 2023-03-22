%% generateEllipsePoints
% Perry Hong 
% 22 March 2023
%
% Generate set of x and y points for plotting of ellipses
% Inputs arguments are all L x 1 column vectors, where L refers to the number of input ellipses
% Angle in degrees, defined ACW from x-axis
% 
% === Optional arguments ===
% res: resolution of parametric theta (default 0.01 radians)

%% Begin function
function [x, y] = generateEllipsePoints(center_x, center_y, major, minor, angle, options)

    arguments
    
        center_x (:,1) double
        center_y (:,1) double
        major (:,1) double
        minor (:,1) double
        angle (:,1) double
        options.res double = 0.01

    end

    angle = deg2rad(angle);

    % Parametric
    theta = 0:options.res:2*pi;

    x = major*cos(theta).*cos(angle) - minor*sin(theta).*sin(angle) + center_x;
    y = major*cos(theta).*sin(angle) + minor*sin(theta).*cos(angle) + center_y;

end

