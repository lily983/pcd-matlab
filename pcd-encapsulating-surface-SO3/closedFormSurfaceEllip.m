%% Loading rotations from samples
clc;clear;close all

load('samples.mat')
R = csv2group(samples, 'SO3');
R = filterData(R, eye(3));
mu = get_mean_cov(R, 'SO', 1);

R_filter = filterData(R, mu);
% R_filter = R;
[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R_filter, 'SO', 1);

%% Loading obj parameters
% axes = [0.034,0.021, 0.016];
axes = rand(1,3)*0.1+0.01;
obj = SuperQuadrics({axes, rand(1,2)*2, [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

%% Get encapsulating surface points
x_new_4 = idea4scale(obj, R_filter, size(R_filter, 3)/1.2);

%%
axes_scale = obj.a * 1.2;
obj_scale = SuperQuadrics({axes_scale, obj.eps, [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});

%% Plot obj at mean rotation and its slightly rotated copies 
% close all;
figure; hold on; axis equal
% obj_mean = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});
% obj_mean.PlotShape('b', 0.1,0.1);

[R_extreme, d_R] = getExtremeR(R_filter, mu_SO3, Sigma_SO3);

for i=1:size(R_extreme,3)/3
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], zeros(3,1), rotm2quat(R_extreme(:,:,i)), [20,20]});
    obj_i.PlotShape('b', 0.001,0.15);
    pause(0.1);
end
% obj_scale.PlotShape(hex2rgb('EBBF00'), 0.3,0.15);% ginger yellow
% plotSurface(x_new_7, obj.N, hex2rgb('EBBF00'), 0.4); % dark green
plotSurface(x_new_4, obj.N, hex2rgb('ff0010'), 0.4); % dark green
%%
figure; hold on; axis equal
% plotSurface(x_new_1, obj.N,  hex2rgb('45AC59')); % dark green
% plotSurface(x_new_2, obj.N, hex2rgb('98646b'));% grape
% plotSurface(x_new_3, obj.N, hex2rgb('2023c7')); % bright purple
% obj.PlotShape('b', 0.5,0.15);
plotSurface(x_new_4, obj.N, hex2rgb('45AC59'), 0.3); %  dark green
% plotSurface(x_new_7, obj.N, hex2rgb('EBBF00'), 0.5); %  grape
plotSurface(x_new_8, obj.N, hex2rgb('ff0010'), 0.3); %  grape
% plotSurface(x_new_8, obj.N, hex2rgb('ff0010')); %  bright purple
% plotSurface(mu_SO3*x_new_51, obj.N, hex2rgb('ff0010')); % red
% plotSurface(mu_SO3*x_new_52, obj.N, hex2rgb('45AC59')); %  dark green
