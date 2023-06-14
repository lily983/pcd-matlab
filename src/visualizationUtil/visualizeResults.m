function visualizeResults(results, objectType)
% visualizeResults Visualize PCD results from each methods
% 
% Inputs
%     results     a struct array containing flag, dist, PCD methods results
%     objectType   if testing objects are sphere, then plot PCD-maxpdf

% create a map container, keySet is the name of PCD-methods, value set is
% the color and data array 
keySet = {'PCD-exact', 'PCD-maxpdf', 'PCD-convex', 'PCD-EB', 'PCD-SG', 'PCD-GMM'};
colorSet = {hex2rgb('45498C')./255, hex2rgb('EE6B31')./255, hex2rgb('45AC59')./255,...
    hex2rgb('EBBF00')./255, hex2rgb('CE2836')./255, hex2rgb('4FAAD1')./255}
dataSet = {results.PCDExactArray, results.PCDMaxpdfArray, results.PCDConvexArray,...
    results.PCDEBArray, results.PCDSGArray, results.PCDGMMArray};
colorMap = containers.Map(keySet, colorSet);
dataMap = containers.Map(keySet, dataSet);

if strcmp(objectType, 'sphere')==0
    keySet(2) = [];
    remove(colorMap, 'PCD-maxpdf');
    remove(dataMap, 'PCD-maxpdf');
end

%% First plot data point and then plot fitting curve
dataNumber = size(results.PCDExactArray,2);
xIndex = 1:1:dataNumber;
markerAlpha = 0.4;
markerSize = 30;

figure; hold on;
ylim([0 1])
xlim([1 dataNumber])

for iMethod=1:dataMap.Count

    currentDataArray = dataMap(keySet{iMethod});
    currentDataColor =  colorMap(keySet{iMethod});
    
    scatter(xIndex, currentDataArray,...
         'MarkerFaceColor', currentDataColor,...
         'AlphaData', currentDataArray,...
         'MarkerFaceAlpha', 'flat',...
         'MarkerEdgeAlpha', 0.0);
    scatter(xIndex, currentDataArray,...
        markerSize,...
         'MarkerFaceColor', currentDataColor,...
         'MarkerFaceAlpha',markerAlpha,...
         'MarkerEdgeAlpha', 0.0);
    visualizeDataFittingResult(currentDataArray, currentDataColor);
end
hold off;
end