%% heatmap2Dpdf_Gaussian
% Perry Hong
% 8 September 2022
%
% Full: Does full gaussian integral for each grid and adds that value to that grid on the heatmap (EXACT)
% Sampled: Sampled gaussian integral 
%
% Output: N x M vector where N and M are the number of bins for Lat and Long. Sum normalised to 1 if desired.
%
% ==== Input arguments ====
% input = L x 5 vector 
% [ellipseLat ellipseLong ellipseMajor ellipseMinor ellipseAngle]
% ellipseMajor and ellipseMinor (95% 2D confidence) in m
% ellipseAngle ACW from x-axis (major is x-axis) in degrees
% nBins = [nBinsLat nBinsLong]
% rangeBins = [rangeLatLower rangeLatUpper; rangeLongLower rangeLongUpper]
%
% === Optional arguments ===
% method = integration method: 'full', 'sampled'
% maxGrid = number of grid points to compute (refer to documentation in code)
% norm = boolean inidicating normalisation of output heatmap

%% Begin function

function heatmap2D_Gaussian = heatmap2Dpdf_Gaussian(input, nBins, rangeBins, options)

    arguments
    
        input (:,5) double
        nBins (1,2) double
        rangeBins (2,2) double
        options.method {mustBeMember(options.method,['full','sampled'])} = 'full'
        options.maxGrid int16 = 15
        options.norm logical = 0 
        
    end

    nBinLat = nBins(1);
    nBinLong = nBins(2);

    rangeBinLat = rangeBins(3,:);
    rangeBinLong = rangeBins(4,:);

    % Check if input is valid (within rangeBins)
    if (~isempty(find(input(:,1)<rangeBinLat(1))) || ~isempty(find(input(:,1)>rangeBinLat(2))))
        error('Input has a value outside of defined latitude range.');
    elseif (~isempty(find(input(:,2)<rangeBinLong(1))) || ~isempty(find(input(:,2)>rangeBinLong(2))))
        error('Input has a value outside of defined longitude range.');
    end

    % Resolution of each bin
    resBinLat = (rangeBinLat(2)-rangeBinLat(1))/nBinLat;
    resBinLong = (rangeBinLong(2)-rangeBinLong(1))/nBinLong;

    L = size(input, 1);

    % Shift latitude and longitude values to [0 x]
    shift = [ones(L, 1)*rangeBinLat(1) ones(L,1)*rangeBinLong(1) zeros(L,1) zeros(L,1) zeros(L,1)];
    inputShifted = input - shift;

    % Initialise output
    heatmap2D_Gaussian = zeros(nBinLat, nBinLong);
    
    % Evaluate each ellipse
    progbar = waitbar(0, 'Generating heatmap...');
    
    for n = 1:L

        tEllipseStart = tic;
        
        currentEllipseLat = inputShifted(n,1);
        currentEllipseLong = inputShifted(n,2);
        currentEllipseMajorRaw = inputShifted(n,3);
        currentEllipseMinorRaw = inputShifted(n,4);
        currentEllipseAngleRaw = inputShifted(n,5);
        
        % Bins are defined w.r.t. to the top right = bottom left bin is (1, 1)
        idxBinLat = ceil(currentEllipseLat/resBinLat);
        idxBinLong = ceil(currentEllipseLong/resBinLong);  

        % Edge cases: should correspond to index 1 or binMax as well
        if idxBinLat == 0 idxBinLat = 1; end
        if idxBinLong == 0 idxBinLong = 1; end
        if idxBinLat == nBinLat+1 idxBinLat = nBinLat; end
        if idxBinLong == nBinLong+1 idxBinLong = nBinLong; end
        
        % For this covariance matrix calculation, angle is defined as the angle the semi-major axis makes with the x-axis, anti-clockwise; semi-major axis is defined as the x-axis
        currentEllipseAngle = deg2rad(currentEllipseAngleRaw);
 
        % Convert ellipseMajor (long) and ellipseMinor (lat) to LL units
        % Note: convertDisttoLL only properly does the conversion along the afixed latitude/longitude - distances at angles dont make sense as 1 degree lat != 1 degree long
        % Here, our semi-major/minor axes are at an angle - we have neglected this and arbitrarily chosen to convert the semi-major along the longitude and semi-minor along the latitude
        % For the error ellipse sizes we are working with (~20km max), the error is small
        trueEllipseLat = currentEllipseLat + rangeBinLat(1);                                
        [currentEllipseMinor, currentEllipseMajor] = convertDisttoLL(currentEllipseMinorRaw, currentEllipseMajorRaw, trueEllipseLat);

        % sigma11 and sigma22 in diagonal basis for the ellipse which is 0.95 cdf
        sigma1diag = currentEllipseMajor/sqrt(chi2inv(0.95, 2));
        sigma2diag = currentEllipseMinor/sqrt(chi2inv(0.95, 2));

        % Compute covariance matrix in LL units -
        % Angle is defined as angle between x-axis and sigma11-axis (anticlockwise)
        sigmaDiag = [sigma1diag^2 0; 0 sigma2diag^2];
        rotMatrix = [cos(currentEllipseAngle) -sin(currentEllipseAngle); sin(currentEllipseAngle) cos(currentEllipseAngle)];
        sigmaRot = (rotMatrix)*sigmaDiag*(rotMatrix.'); % sigmaRot = R sigmaDiag R_T 
      
        mu = [currentEllipseLong currentEllipseLat];

        %% Evaluation range - dynamic evaluation range based on ellipse
        % semi-major axis relative to the grid size, we at least compute
        % until 95% of the Gaussian (which conveniently defines the ellipse in this
        % case, so we just use the ellipse semi-major value), or maximum of
        % maxGrid grids away (compute a maximum 11 x 11 matrix if maxGrid = 5)
        % adjust accordingly based on gridsize vs. expected semi-major 
        
        checkGridRange = min([options.maxGrid max(ceil([currentEllipseMajor/resBinLat currentEllipseMajor/resBinLong]))]); % accounts for arbitrary rotation of ellipse by just using semi-major axis 
        gridSize = checkGridRange*2+1;
        gridCDF = zeros(gridSize,gridSize);

        % Full integration over entire grid
        if strcmp(options.method, 'full')
            
            for idxLat = 1:gridSize
                
                for idxLong = 1:gridSize

                    idxCurrentGridLat = idxLat - checkGridRange - 1;        % which grid we are looking at now (-3 -2 -1 0 1 2 3)
                    idxCurrentGridLong = idxLong - checkGridRange - 1;

                    currentGridLat = (idxBinLat+idxCurrentGridLat)*resBinLat;             % LL value of the top right of the current grid
                    currentGridLong = (idxBinLong+idxCurrentGridLong)*resBinLong;
                    
                    gridCDF(idxLat, idxLong) = max([0 mvncdf([currentGridLong-resBinLong currentGridLat-resBinLat], [currentGridLong currentGridLat], mu, sigmaRot)]);  % gridSize-idxLat+1 because array rows indexes from top to bottom
                   
                    %gridSize-idxLat+1
                    
                end % end of loop over grid longitudes

            end % end of loop over grid latitudes
 
        % Sampled integration over subpoints on grid
        elseif strcmp(options.method, 'sampled')
        
            subGridNum = 3;
            subGridRes = 1/(subGridNum+1); 
        
            % Generate CDF grid for this ellipse 
            for idxLat = 1:gridSize

                for idxLong = 1:gridSize

                    idxCurrentGridLat = idxLat - checkGridRange - 1;        % which grid we are looking at now (-3 -2 -1 0 1 2 3)
                    idxCurrentGridLong = idxLong - checkGridRange - 1;
                    currentGridSum = 0;

                    % Evaluate PDF at each subpoint on the grid
                    for idxSubGridLat = 1:subGridNum

                        for idxSubGridLong = 1:subGridNum

                            currentGridPointLat = (idxBinLat+idxCurrentGridLat-(idxSubGridLat*subGridRes))*resBinLat;
                            currentGridPointLong = (idxBinLong+idxCurrentGridLong-(idxSubGridLat*subGridRes))*resBinLong;

                            currentGridSum = currentGridSum + max([0 mvnpdf([currentGridPointLong currentGridPointLat], mu, sigmaRot)]);
                            
                        end % end of loop over subgrid longitudes

                    end % end of loop over subgrid latitudes

                    gridCDF(idxLat, idxLong) = currentGridSum/(subGridNum^2);   % gridSize-idxLat+1 because array rows indexes from top to bottom
                    
                end % end of loop over grid longitudes

            end % end of loop over grid latitudes

           gridCDF = gridCDF/sum(sum(gridCDF));    % normalisation required, as mvnpdf is not normalised to 1
            
        end % end of CDF calculation
         
        %% Delete rows and columns which exceed boundaries of LL bins
        excessLeft = idxBinLong - (checkGridRange + 1);         % extra +1 for the left and bottom cases is due to the convention of indexing of the grids, starts at (0 0)
        excessRight = nBinLong - (idxBinLong + checkGridRange);
        excessTop = nBinLat - (idxBinLat + checkGridRange);
        excessBottom = idxBinLat - (checkGridRange + 1);
        
        nLeft = checkGridRange; nRight = checkGridRange; nTop = checkGridRange; nBottom = checkGridRange;
        
        if excessLeft < 0       
            gridCDF = gridCDF(:,1+abs(excessLeft):end);   % delete left abs(excessLeft) columns
            nLeft = nLeft + excessLeft; 
        end          
        if excessRight < 0      
            gridCDF = gridCDF(:,1:end-abs(excessRight));  % delete right abs(excessRight) columns
            nRight = nRight + excessRight; 
        end        
        if excessTop < 0        
            gridCDF = gridCDF(1:end-abs(excessTop),:);    % delete left abs(excessTop) rows
            nTop = nTop + excessTop; 
        end           
        if excessBottom < 0     
            gridCDF = gridCDF(1+abs(excessBottom):end,:); % delete right abs(excessBottom) rows
            nBottom = nBottom + excessBottom; 
        end        
        
        % Add grid values to heatmap
        heatmap2D_Gaussian(idxBinLat-nBottom:idxBinLat+nTop, idxBinLong-nLeft:idxBinLong+nRight) = heatmap2D_Gaussian(idxBinLat-nBottom:idxBinLat+nTop, idxBinLong-nLeft:idxBinLong+nRight) + gridCDF;
        
        waitbar(n/L, progbar, ['Generating heatmap: computing ellipse ' num2str(n) ' of ' num2str(L)]);
        disp(['Ellipse n=' num2str(n) ' took ' num2str(toc(tEllipseStart)) 's to evaluate.']);

    end % end of this current ellipse data point
    
    close(progbar);
    
    % Normalise if desired
    if options.norm
        heatmap2D_Gaussian = heatmap2D_Gaussian/(sum(sum(heatmap2D_Gaussian));
    end
    
end


