clc;close all;clear

%% Loading rotations from samples
R = csv2group(samples, 'SO3');

R_filter = filterData(R, mu_SO3);

%% Generate ellipsoid
testNum = 100;

axesArray = 0.01 + (0.1 - 0.01) * rand(testNum, 3);

axesArray = sort(axesArray, 2, 'descend');

d_5_1 = zeros(1, testNum);
d_5_3 = zeros(1, testNum);

%%
for i=1:testNum
    obj_i = SuperQuadrics({axesArray(i, :), [1,1], [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});
    
    [~, d_5_1(i)] = idea5(obj_i, R, '1' );
    
    [~, d_5_3(i)] = idea5(obj_i, R, '3' );
    
end
%%
difference_1_3 = d_5_1- d_5_3;
[max_diff, max_index] = max(difference_1_3)

[min_diff, min_index] = min(difference_1_3)

[~, mean_index] = min(abs(difference_1_3 -  mean(difference_1_3)))
mean_diff  = difference_1_3(mean_index)

[~, median_index] = min(abs(difference_1_3 -  median(difference_1_3)))
median_diff  = difference_1_3(median_index)

%%
figure; hold on; axis equal
[max_obj, max_obj_maxd]=visualizeEncapsulatingSurface(axesArray, max_index, '3', R_filter, mu_SO3,  1,  hex2rgb('EBBF00'),0.8);
visualizeEncapsulatingSurface(axesArray, max_index, '3', R_filter, mu_SO3, 1,  hex2rgb('45AC59'),0.2);
max_obj.PlotShape('b', 0.3);
max_obj_maxd.PlotShape('r', 0.1);
% title('Original ellip (blue), surface of improvement #2 (green), and surface of improvement #1 (yellow)')

%%
figure; hold on; axis equal
[min_obj, min_obj_maxd]=visualizeEncapsulatingSurface(axesArray, min_index, '1', R_filter, mu_SO3, 1,  hex2rgb('EBBF00'),0.2);
visualizeEncapsulatingSurface(axesArray, min_index, '3', R_filter, mu_SO3, 1,  hex2rgb('45AC59'),0.3);
min_obj.PlotShape('b', 0.3);
min_obj_maxd.PlotShape('r', 0.1);

% title('Original ellip (blue), surface of improvement #2 (green), and surface of improvement #1 (yellow)')
%%
figure; hold on; axis equal
[mean_obj, mean_obj_maxd1, ~, mean_d_1]=visualizeEncapsulatingSurface(axesArray, mean_index, '3', R_filter, mu_SO3, 1,  hex2rgb('EBBF00'),0.8);
% [median_obj, median_obj_maxd3, ~, median_d_3]=visualizeEncapsulatingSurface(axesArray, mean_index, '3', R_filter, mu_SO3, 1,  hex2rgb('45AC59'),0.8);
mean_obj.PlotShape('b', 0.3);
mean_obj_maxd1.PlotShape('r', 0.1);
% title('Original ellip (blue), surface of improvement #2 (green), and surface of improvement #1 (yellow)')
%%
figure; hold on; axis equal
[median_obj, median_obj_maxd_1, ~, median_d_1]=visualizeEncapsulatingSurface(axesArray, median_index, '1', R_filter, mu_SO3, 1,  hex2rgb('EBBF00'),0.8);
% [~, median_obj_maxd_3, ~, median_d_3]=visualizeEncapsulatingSurface(axesArray, median_index, '3', R_filter, mu_SO3,   1,  hex2rgb('45AC59'),0.3);
median_obj.PlotShape('b', 0.3);
median_obj_maxd_1.PlotShape('b', 0.1);
% title('Original ellip (blue), surface of improvement #2 (green), and surface of improvement #1 (yellow)')

%%
figure; hold on 
% plot(d_5_1, 'g')
scatter(1:testNum, d_5_1, 'b')
scatter(1:testNum, d_5_3, 'r')

%%
figure; hold on 
boxplot(difference_1_3')
ylim([-0.001, 0.03])
title('Difference between improvement 2 and improvement 1')

%%
function [obj_mu, obj_maxd,x_new, d]=visualizeEncapsulatingSurface(axesarray, index, method, R, mu_SO3, isplot, color, alpha)

if nargin == 5
    isplot=false;
end

axes = axesarray(index,:);

obj = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

obj_mu = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});

[x_new, d, d_index] = idea5(obj, R, method );

obj_maxd = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,d_index)), [20,20]});

if isplot
    plotSurface(mu_SO3 * x_new, [20,20], color, alpha); % red
    
    % plot a*R*e1
    e1 = [1;0;0];
    e1 = axes(1) * R(:,:,d_index) * e1;
    scatter3(e1(1), e1(2), e1(3), 20, 'r');
end

end
