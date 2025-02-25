
clc;clear;close all

%% Loading rotations from samples
T = csv2group(summaryposedata, 'PCG3');
T_filter = filterData(T, eye(3));

R = T_filter(1:3,1:3,:);
mu = get_mean_cov(R, 'SO', 1);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 1);

% R = csv2group(samples, 'SO3');

%  R_filter = filterData(R, mu_SO3);
R_filter = R;

[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R_filter, 'SO', 1);

%% Loading obj parameters
% axes = [0.034,0.021, 0.016];
axes=[0.2927352  0.2915878  0.12707044]';
axes = axes/2
%% wooden bowl
axes = [0.2927352017397958,0.29158779787962896,0.12707043772973667]';
axes = axes/2
%% bottel
axes = [0.26966692120426106,0.18009964067697282,0.10951677121080064]';
axes = axes/2
%% tiny chair
axes = [0.27727629926627784,0.2358111973954244,0.17633288963473515]';
axes = axes/2
%%
% obj = SuperQuadrics({axes, [0.2,0.2], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});
obj = SuperQuadrics({axes, rand(1,2)*2, [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

figure
obj.PlotShape('b', 0.1,0.1);
%% Get encapsulating surface points
x_new_4 = idea4scale(obj, R_filter, size(R_filter, 3)/1.3);

%%
x_new_8 = idea8(obj, R, size(R, 3)/1);

%% Plot obj at mean rotation and its slightly rotated copies 
% close all;
figure; hold on; axis equal

for i=1:size(R_filter,3)
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], zeros(3,1), rotm2quat(R_filter(:,:,i)), [20,20]});
    obj_i.PlotShape('b', 0.001,0.15);
    pause(0.1);
end

plotSurface(x_new_4, obj.N, hex2rgb('45AC59'), 0.3);

%% Plot obj at mean rotation and its slightly rotated copies 
% close all;
figure; hold on; axis equal

for i=1:size(R_filter,3)
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], zeros(3,1), rotm2quat(R_filter(:,:,i)), [20,20]});
    obj_i.PlotShape('b', 0.001,0.15);
    pause(0.1);
end
plotSurface(x_new_8, obj.N, hex2rgb('45AC59'), 0.3); % red
