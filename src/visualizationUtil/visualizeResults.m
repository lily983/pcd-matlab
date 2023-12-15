function visualizeResults(results, keySet, colorSet, sampleNumber)
% visualizeResults Visualize PCD results from each methods

% First plot data point and then plot fitting curve
xIndex = 1:1:sampleNumber;
markerAlpha = 0.5;
markerSize = 30;

figure; hold on;
ylim([0 1])
xlim([1 sampleNumber])

for i = 1:size(keySet,2)
    iMethod = string(keySet(i));
    
    currentDataArray = results.(genvarname(iMethod));
    currentDataColor =  colorSet{i};
    
    scatter(xIndex, currentDataArray,...
         'MarkerFaceColor', currentDataColor,...
         'AlphaData', real(currentDataArray),...
         'MarkerFaceAlpha', 'flat',...
         'MarkerEdgeAlpha', 0.0);
    scatter(xIndex, currentDataArray,...
        markerSize,...
         'MarkerFaceColor', currentDataColor,...
         'MarkerFaceAlpha',markerAlpha,...
         'MarkerEdgeAlpha', 0.0);
     
     if strcmp(iMethod, 'EB_99') || strcmp(iMethod, 'EB_95') || strcmp(iMethod, 'Maxpdf_SQ')
         continue;
     end
    visualizeDataFittingResult(currentDataArray, currentDataColor);
end
plot(xIndex,0.05*ones(sampleNumber),'--');
hold off;
end