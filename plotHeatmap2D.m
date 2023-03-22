%% plotHeatmap2D
% Perry Hong
% 22 March 2023
%
% Plots heatmap on the current figure
% Heatmap generated from heatmap2Dpdf_Gaussian
%
% === Input arguments ===
% heatmap2D_Gaussian = N x M heatmap where N and M are the number of bins for Lat and Long.
% rangeLong: [rangeLongLower rangeLongUpper]
% rangeLat: [rangeLatLower rangeLatUpper]


%% Begin function

function heatmap2D_handle = plotHeatmap2D(heatmap2D_Gaussian, rangeLong, rangeLat)

    arguments
    
        heatmap2D_Gaussian (:,:) double
        rangeLong (1,2) double
        rangeLat (1,2) double
        
    end

    heatmap2D_handle = imagesc('XData',rangeLong,'YData',rangeLat,'CData', heatmap2D_Gaussian); 
    xlim(rangeLong); ylim(rangeLat); 
    colorbar; colormap(jet);

end