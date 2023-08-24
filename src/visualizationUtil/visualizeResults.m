function visualizeResults(results, objectType)
% visualizeResults Visualize PCD results from each methods
% 
% Inputs
%     results     a struct array containing flag, dist, PCD methods results
%     objectType   if testing objects are sphere, then plot PCD-maxpdf

% create a map container, keySet is the name of PCD-methods, value set is
% the color and data array 
keySet = {'PCD-exact', 'PCD-maxpdf', 'PCD-convex', ...
    'PCD-EB-99', 'PCD-GMM-2SG', 'PCD-GMM-4SG', ...
    'PCD-GMM-5SG', 'PCD-EB-95', 'PCD-ellip-exact',...
    'PCD-ellip-bound'};

colorSet = {hex2rgb('45498C'), hex2rgb('EE6B31'), hex2rgb('45AC59'),...
    hex2rgb('EBBF00'), hex2rgb('4FAAD1'), hex2rgb('CE2836'), ...
    hex2rgb('9f45b0'), hex2rgb('98646b'), hex2rgb('b5a642'), ...
    hex2rgb('E6A9CC')};

dataSet = {results.PCDExactArray, results.PCDMaxpdfArray, results.PCDConvexArray,...
    results.PCDEB99Array, results.PCDGMM2SGArray, results.PCDGMM4SGArray, ...
    results.PCDGMM5SGArray, results.PCDEB95Array, results.PCDEllipExactArray,...
    results.PCDEllipBoundArray};

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
markerAlpha = 0.5;
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
     
     if strcmp(keySet{iMethod}, 'PCD-EB-99') || strcmp(keySet{iMethod}, 'PCD-EB-95')...
             || strcmp(keySet{iMethod}, 'PCD-ellip-exact')
         continue;
     end
    visualizeDataFittingResult(currentDataArray, currentDataColor);
end
plot(xIndex,0.05*ones(dataNumber),'--');
hold off;
end