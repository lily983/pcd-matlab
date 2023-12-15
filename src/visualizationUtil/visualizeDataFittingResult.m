function visualizeDataFittingResult(data, color)
% visualizeDataFittingResult Fit data using SLM - Shape Language Modeling and plot fit curve in color
% 
% Inputs
%     data    1xN array of PCD results by a method
%     color   1x3 rgb value
 
% Check NaN
if any(isnan(data))~=0
    return
end

dataSize = size(data,2);

xIndex = 1:1:dataSize;

slm = slmengine(xIndex, data,'plot','off','knots',dataSize,'increasing','on', ...
'leftslope',0,'rightslope',0);

fittingResult = slmeval(xIndex, slm);

plot(xIndex, fittingResult, 'Color', color, 'LineWidth', 2);

end