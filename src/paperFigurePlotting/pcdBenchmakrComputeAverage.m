function [ratio, time]=pcdBenchmakrComputeAverage(results, objectType)
% 
% Inputs
%     results     a struct array containing flag, dist, PCD methods results
%     objectType   if testing objects are sphere, then plot PCD-maxpdf
% create a map container, keySet is the name of PCD-methods, value set is
% the color and data array 
keySet = {'PCD-exact', 'PCD-EB-99', 'PCD-EB-95', ...
    'PCD-ellip-bound', 'PCD-maxpdf', 'PCD-convex', ...
    'PCD-GMM'};

dataSet = {results.PCDExactArray, results.PCDEB99Array, results.PCDEB95Array,...
    results.PCDEllipBoundArray, results.PCDMaxpdfArray, results.PCDConvexArray, ...
    results.PCDGMM2SGArray};

timeSet = {results.exactT, results.EBT, results.EBT,...
    results.EllipBoundT, results.maxpdfT, results.convexT, ...
    results.GMM2SGT};

dataMap = containers.Map(keySet, dataSet);
timeMap = containers.Map(keySet, timeSet);

if strcmp(objectType, 'sphere')==0
    keySet(2) = [];
    remove(dataMap, 'PCD-maxpdf');
end

%% First plot data point and then plot fitting curve
dataNumber = size(results.PCDExactArray,2);
xIndex = 1:1:dataNumber;
markerAlpha = 0.5;
markerSize = 30;

% figure; hold on;
ylim([0 1])
xlim([1 dataNumber])

groundTruthArray=dataMap('PCD-exact');
groungTruthNumber=sum(groundTruthArray<=0.05);

for iMethod=1:dataMap.Count

    currentDataArray = dataMap(keySet{iMethod});
    time.keySet{iMethod} = mean(timeMap(keySet{iMethod}));
    ratio.keySet{iMethod} = sum(currentDataArray<=0.05)/groungTruthNumber;
end


end