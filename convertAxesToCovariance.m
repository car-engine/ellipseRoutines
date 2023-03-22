%% convertAxesToCovariance
% Perry Hong 
% 22 March 2023
%
% Coverts a given ellipse defined by its semi-major, semi-minor and angle
% (ACW w.r.t x-axis, in degrees), return the covariance matrix of its
% underlying Gaussian distribution

%% Begin function
function sigmaRot = convertAxesToCovariance(major, minor, angle, options)

    arguments
    
        major double
        minor double
        angle double
        options.confidence double {mustBeLessThan(options.confidence, 1), mustBeGreaterThan(options.confidence, 0)} = 0.95

    end

    major_scaled = major/sqrt(chi2inv(options.confidence, 2));
    minor_scaled = minor/sqrt(chi2inv(options.confidence, 2));

    t = deg2rad(angle);
    sigmaDiag = [major_scaled^2 0; 0 minor_scaled^2];
    rotMatrix = [cos(t) -sin(t); sin(t) cos(t)];
    sigmaRot = (rotMatrix)*sigmaDiag*rotMatrix.'; % sigmaRot = R sigmaDiag R_T

end
