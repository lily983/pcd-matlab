close all;

% Load the data
load('SingleError_SQ_100_10k.mat');

% Extract the data arrays
baseline = resultsSQ.Exact(:);
lcc_center_point = resultsSQ.LCC_center_point(:);
lcc_tangent_point_cfc = resultsSQ.LCC_tangent_point_cfc(:);
divergence_mesh = resultsSQ.Divergence_mesh(:);

% Calculate the differences
diff_lcc_center_point = lcc_center_point - baseline;
diff_lcc_tangent_point_cfc = lcc_tangent_point_cfc - baseline;
diff_divergence_mesh = divergence_mesh - baseline;

% Prepare the data for the violin plot
% data = {diff_lcc_center_point, diff_lcc_tangent_point_cfc, diff_divergence_mesh}';
% labels = {'LCC Center Point', 'LCC Tangent Point CFC', 'Divergence Mesh'};
% labels = cellstr(labels);
data = [diff_lcc_center_point, diff_lcc_center_point];
labels = {'LCC Center Point', 'LCC Tangent Point CFC'};

% Create the violin plot
figure; 
hold on

vs = violinplot(data, labels);

vs(1,1).ShowBox = 1;
vs(1,1).QuartileStyle = 'shadow';
vs(1,1).HalfViolin = 'left';
vs(1,1).ShowMean =0;
% % vs(1).MarkerSize = 5;


% Customize the violin plot
% for i = 1:length(vs)
%     % Set the face color of the violin
% %     vs(i).ViolinColor = [0.5 0.5 0.5]; % Gray color for the violin body
% 
%     % Set the color for the mean line
%     vs(i).MeanPlot.Color = 'b'; % Blue color for the mean
%     vs(i).MeanPlot.LineWidth = 8;
% 
%     % Set the color for the median line
%     vs(i).MedianPlot.Color = 'r'; % Red color for the median
%     vs(i).MedianPlot.LineWidth = 8;
%     
%     % Set the color for the quartile lines
%     vs(i).BoxColor = 'g'; % Green color for the quartiles
%     vs(i).BoxWidth = 8;
% end

title('Violin Plot of Differences from Exact with Quartile Lines');
ylabel('Difference');
grid on;
