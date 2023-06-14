function visualizeTwoMethodsIntersection(firstIndex, secondIndex, PCDExactDataArray)
% visualizeTwoMethodsIntersection Plot the intersection point of two PCD
% methods results array. In R3 case, call this function to emphasize the
% diffenrece between PCD-convex and PCD-GMM
% 
% Inputs
%     firstIndex  index of the first intersection point
%     secondIndex  index of the second intersection point
%     PCDExactDataArray  PCD-exact array, which is specified on y-axis

prob=PCDExactDataArray(firstIndex);
prob2=PCDExactDataArray(secondIndex);

x1_points = [0, 0, firstIndex, firstIndex];  
y1_points = [0, prob, prob, 0];
color1 = hex2rgb('06d6a0')./255;
a1 = fill(x1_points, y1_points, color1);
a1.FaceAlpha = 0.6;
a1.EdgeAlpha = 0.0;

x2_points = [firstIndex, firstIndex, secondIndex, secondIndex];  
y2_points = [prob, 1, 1, prob];
color2 =  hex2rgb('118ab2 ')./255;
a2 = fill(x2_points, y2_points, color2);
a2.FaceAlpha = 0.2;
a2.EdgeAlpha = 0.0;
% 
vertical_line = linspace(0, 1, 100);
plot(firstIndex*ones(size(vertical_line)), vertical_line, 'Color', hex2rgb('073b4c')./255);

horizontal_line = prob * ones(1, secondIndex);
plot(1:secondIndex, horizontal_line, 'Color', hex2rgb('073b4c')./255);

vertical_line2 = linspace(prob, 1, 100);
plot(secondIndex*ones(size(vertical_line)), vertical_line2, 'Color', hex2rgb('073b4c')./255);

horizontal_line2 = prob2 * ones(1, secondIndex);
plot(1:secondIndex, horizontal_line2, 'Color', hex2rgb('073b4c')./255, 'LineStyle',":", 'LineWidth', 2);

horizontal_line3 = prob * ones(1, N-secondIndex+1);
plot(secondIndex:N, horizontal_line3, 'Color', hex2rgb('073b4c')./255, 'LineStyle',":", 'LineWidth', 2);

horizontal_line4 = prob2 * ones(1, N-secondIndex+1);
plot(secondIndex:N, horizontal_line4, 'Color', hex2rgb('073b4c')./255, 'LineStyle',":", 'LineWidth', 2);

yticks([0 prob  0.2  0.4 prob2  0.6  0.8  1.0])
ylim([0 1])
xlim([1 500])

end