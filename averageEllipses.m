%% averageEllipses
% Perry Hong
% 8 Sep 2022
%
% Computes averaging of ellipses on WGS84 using two methods
% Davis: Combining Error Ellipses, Davis J.E. (2007) - cxc.cda.harvard.edu/csc/memos/files/Davis_ellipse.pdf
% Berkeley: Data Analysis Toolkit 12: Weighted averages and their uncertainties, Kirchner J (2006)
% semi-major = x-axis
% angle is defined as the angle semi-major axis makes wrt to x-axis, anti-clockwise
% Inputs: mu = [long, lat]; major, minor (95% 2D confidence) in meters; angle in degrees

%% Begin function
function [mu_weighted_mean, sigma_davis, sigma_berkeley] = averageEllipses(mu, major, minor, angles)

    arguments

        mu (:,2) double
        major (1,:) double
        minor (1,:) double
        angles (1,:) double

    end

    L = length(major);
    minor_LL = zeros(L,1);
    major_LL = zeros(L,1);

    % Convert to covariance matrix 1 sigma
    major = major/sqrt(chi2inv(0.95, 2));
    minor = minor/sqrt(chi2inv(0.95, 2));

    % Convert ellipseMajor (long) and ellipseMinor (lat) to LL units
    % Note: convertDisttoLL only properly does the conversion along a fixed latitude/longitude - distances at angles dont make sense as 1 degree lat != 1 degree long
    % Here, our semi-major/minor axes are at an angle - we have neglected this and arbitrarily chosen to convert the semi-major along the longitude and semi-minor along the latitude
    % For the error ellipse sizes we are working with, the error is small
    for idx = 1:L 
        [minor_LL(idx,1), major_LL(idx,1)] = convertDisttoLL(minor(idx), major(idx), mu(idx,2));
    end
   
    % Calculate covariance matrices
    sigmaRot = cell(1,L);
    for n = 1:L
        t = deg2rad(angles(n));
        sigmaDiagCurrent = [major_LL(n)^2 0; 0 minor_LL(n)^2]; % rotated frame diagonal matrix
        rotMatrix = [cos(t) -sin(t); sin(t) cos(t)];
        sigmaRot{n} = (rotMatrix)*sigmaDiagCurrent*rotMatrix.'; % sigmaRot = R sigmaDiag R_T
    end

    %% MIT Davis Ellipse
    sigmaRotInvSum = zeros(2,2);

    for n = 1:L 
        sigmaRotInvSum = sigmaRotInvSum + inv(sigmaRot{n});
    end

    sigma_davis = inv(sigmaRotInvSum);

    % Calculate muAvg 
    sumInvCovMu = zeros(2,1);

    for n = 1:L
       sumInvCovMu = sumInvCovMu + sigmaRot{n}\(mu(n,:).');
    end

    mu_weighted_mean = sigma_davis*sumInvCovMu;
    
    %% Berkeley Method
    
    weights = cell(1,L);
    
    for n = 1:L
       weights{n} = sigma_davis/sigmaRot{n};
    end
    
    numerator = zeros(2,2);

    for n = 1:L
        numerator = numerator + weights{n}*((mu(n,:).'-mu_weighted_mean)*(mu(n,:).'-mu_weighted_mean).');  
    end
    
    sigma_berkeley = numerator*size(mu,1)/(size(mu,1)-1)/size(mu,1);
    
end
