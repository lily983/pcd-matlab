
clc;close all

%% Loading rotations from samples
T = csv2group(marker38new, 'PCG3');
% [T_filter, filter_indices] = filterData(T, eye(3));
T_filter = T;

R = T_filter(1:3,1:3,:);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 1);


[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R, 'SO', 1);

[T_filter, filtered_indices] = filterData(T, eye(3));

% T_filter = T;

R = T_filter(1:3,1:3,:);
tc = T_filter(1:3,4,:);
mu = get_mean_cov(R, 'SO', 0);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 0);

[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R, 'SO', 0);

[mu_t(1:3,4)', rotm2axang(mu_SO3)]

%% wooden bowl
axes = [0.2927352017397958,0.29158779787962896,0.12707043772973667]';
axes = axes/2
%% tiny chair
axes = [0.27727629926627784,0.2358111973954244,0.17633288963473515]';
axes = axes/2
%% container box
axes=[0.06,0.3,0.065]';
%% form board
axes=[0.205,0.28,0.0225]';
%%
% obj = SuperQuadrics({axes, [0.2,0.2], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});
obj = SuperQuadrics({axes, [0.2,0.2], [0,0], mu_t(1:3,4), rotm2quat(mu_SO3), [20,20]});

figure; axis equal
obj.PlotShape('b', 0.1,0.1);
%% Get encapsulating surface points
x_new_4 = idea4scale(obj, R_filter, size(R_filter, 3)/1.1);

%%

x_new_8 = idea8(obj, R_filter, size(R_filter, 3)/1.1);

%%
axes_scale = obj.a * 1.3
obj_scale = SuperQuadrics({axes_scale, obj.eps, [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});

%%
scale = 1.4;

sub_obj = EnlargedSuperQuadrics(obj, R, scale);

%%
figure; hold on; axis equal
obj.PlotShape('g', 0.3);
for i=1:size(R,3)
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], tc(1:3,i), rotm2quat(R(:,:,i)), [20,20]});
    obj_i.PlotShape('g', 0.05,0.05);
    pause(0.1);
end
sub_obj.PlotShape('b', 0.2);

%% Plot obj at mean rotation and its slightly rotated copies 
% close all;
figure; hold on; axis equal

for i=1:size(R_filter,3)/2
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], zeros(3,1), rotm2quat(R_filter(:,:,i)), [20,20]});
    obj_i.PlotShape('b', 0.001,0.15);
    pause(0.1);
end
% % 
obj_scale.PlotShape(hex2rgb('EBBF00'), 0.2,0.1);% ginger yellow
plotSurface(x_new_8, obj.N, hex2rgb('45AC59'), 0.3); % red
% plotSurface(x_new_4, obj.N, hex2rgb('45AC59'), 0.3); % green
%%
figure; hold on; axis equal;
% obj_scale.PlotShape(hex2rgb('EBBF00'), 0.2,0.15);
plotSurface(x_new_4, obj.N, hex2rgb('45AC59'), 0.1);
plotSurface(x_new_8, obj.N, hex2rgb('98646b'), 0.1); % red

%% 

test_esq = EnlargedSuperQuadrics(obj, R, 1.1);
