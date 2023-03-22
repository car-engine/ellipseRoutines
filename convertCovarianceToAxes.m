%% convertCovarianceToAxes
% Perry Hong 
% 22 March 2023
%
% Given the covariance matrix of a Gaussian distribution, output the
% semi-major, semi-minor, and angle (ACQ w.r.t x-axis, in degrees) of an ellipse satisfying 
% a given confidence interval

%% Begin function
function [major, minor, angle] = convertCovarianceToAxes(sigmaRot, options)

    arguments

        sigmaRot (2,2) double
        options.confidence double {mustBeLessThan(options.confidence, 1), mustBeGreaterThan(options.confidence, 0)} = 0.95

    end

    [rot, ellipse_eig, ~] = svd(sigmaRot);

    major = sqrt(ellipse_eig(1,1))*sqrt(chi2inv(options.confidence,2));
    minor = sqrt(ellipse_eig(2,2))*sqrt(chi2inv(options.confidence,2));
    angle = rad2deg(atan2(rot(2,1),rot(1,1)));

end