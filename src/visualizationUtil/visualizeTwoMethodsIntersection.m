function visualizeTwoMethodsIntersection(firstIndex, secondIndex, PCDExactArray)
% visualizeTwoMethodsIntersection Plot the intersection point of two PCD
% methods results array. In R3 case, call this function to emphasize the
% diffenrece between PCD-convex and PCD-GMM
% 
% Inputs
%     firstIndex  index of the first intersection point
%     secondIndex  index of the second intersection point
%     PCDExactArray  PCD-exact array, which is specified on y-axis
hold on;
prob=PCDExactArray(firstIndex);
prob2=PCDExactArray(secondIndex);

firstPoint = [0, 0, firstIndex, firstIndex];  
firstPointProb = [0, prob, prob, 0];
color1 = hex2rgb('06d6a0');
area1 = fill(firstPoint, firstPointProb, color1);
area1.FaceAlpha = 0.5;
area1.EdgeAlpha = 0.0;

secndPoint = [firstIndex, firstIndex, secondIndex, secondIndex];  
secondPointProb = [prob, 1, 1, prob];
color2 =  hex2rgb('118ab2 ');
area2 = fill(secndPoint, secondPointProb, color2);
area2.FaceAlpha = 0.2;
area2.EdgeAlpha = 0.0;

arraySize = max(size(PCDExactArray, 2), size(PCDExactArray, 2));
verticalLine1 = linspace(0, 1, 100);
plot(firstIndex*ones(size(verticalLine1)), verticalLine1, 'Color', hex2rgb('073b4c'));

horizontalLine1 = prob * ones(1, secondIndex);
plot(1:secondIndex, horizontalLine1, 'Color', hex2rgb('073b4c'));

verticalLine2 = linspace(prob, 1, 100);
plot(secondIndex*ones(size(verticalLine1)), verticalLine2, 'Color',  hex2rgb('073b4c'));

horizontalLine2 = prob2 * ones(1, secondIndex);
plot(1:secondIndex, horizontalLine2, 'Color', hex2rgb('073b4c'), 'LineStyle',":", 'LineWidth', 2);

hotizontalLine3 = prob * ones(1, arraySize-secondIndex+1);
plot(secondIndex:arraySize, hotizontalLine3, 'Color', hex2rgb('073b4c'), 'LineStyle',":", 'LineWidth', 2);

horizontalLine4 = prob2 * ones(1, arraySize-secondIndex+1);
plot(secondIndex:arraySize, horizontalLine4, 'Color', hex2rgb('073b4c'), 'LineStyle',":", 'LineWidth', 2);

yticks(sort([0 0.2 0.4 0.6 0.8 1.0 prob  prob2], 2))
ylim([0 1])
xlim([1 arraySize])

end