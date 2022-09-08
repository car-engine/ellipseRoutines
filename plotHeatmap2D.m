%% plotHeatmap2D
% Plots heatmap on the current figure
% Heatmap generated from heatmap2Dpdf_Gaussian
%
% === Input arguments ===
% heatmap2D_Gaussian = N x M heatmap where N and M are the number of bins for Lat and Long
% rangeLat: [rangeLatLower rangeLatUpper]
% rangeLong: [rangeLongLower rangeLongUpper]

%% Begin function

function heatmap_handle = plotHeatmap2D(heatmap2D_Gaussian, rangeLat, rangeLong)

    heatmap_handle = imagesc('XData',rangeLong,'YData',rangeLat,'CData', heatmap2D_Gaussian); 
    ylim(rangeLat); xlim(rangeLong);
    colorbar; colormap(jet);

end