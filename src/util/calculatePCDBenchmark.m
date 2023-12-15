function [ratio, time, upperBound] = calculatePCDBenchmark(results, keySet, TwoErrorCase)
%% Calculate the results of benchmark and return two metrics. ratio: ratio of available points, time: mean of computation time
if TwoErrorCase
    groundTruthArray  = results.Exact_two;
else 
    groundTruthArray  = results.Exact;
end

% Get number of points below threshold (0.05) of the ground truth array
groungTruthNumber=sum(groundTruthArray<=0.05);

for i = 1:size(keySet, 2)
    iMethod = string(keySet(i));
    currentDataArray = results.(genvarname(iMethod));
    currentTimeArray = results.(genvarname(append(iMethod+'Time')));
    % ratio of points inside threshold
    ratio.(genvarname(iMethod)) = sum(currentDataArray<=0.05)/groungTruthNumber;
    % ratio of upper bound
    upperBound.(genvarname(iMethod)) = sum(currentDataArray - groundTruthArray < - 0.001)/size(groundTruthArray, 2);
%     remove NaN in time array
    if all(isnan(currentTimeArray))
        time.(genvarname(iMethod)) = NaN;
    else
        time.(genvarname(iMethod)) = mean(rmmissing(currentTimeArray));
    end
end

end