
clc;close all

%% Loading rotations from samples
T = csv2group(tinychair, 'PCG3');
% T = csv2group(woodenbowl1, 'PCG3');
[T_filter, filter_indices] = filterData(T, eye(3));
% T_filter = T;
get_mean_cov(T_filter, 'PCG', 1);

%  print('s1s2positionerror', '-dpng', '-r600')

R = T_filter(1:3,1:3,:);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 1);

[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R, 'SO', 1);

%%
for i=1:size(R,3)
    t(i,:) = T_filter(1:3,4,i);
    quat(i,:) = rotm2quat(R(:,:,i));
    pose_csv(i,:) = [t(i,:), quat(i,2:4), quat(i,1)];
end

%% wooden bowl
axes = [0.2927352017397958,0.29158779787962896,0.12707043772973667]';
axes = axes/2

%% tiny chair
axes = [0.27727629926627784,0.2358111973954244,0.17633288963473515]';
axes = axes/2
%%
% obj = SuperQuadrics({axes, [0.2,0.2], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});
obj = SuperQuadrics({axes, [0.2,0.2], [0,0], mu_t(1:3,4), rotm2quat(mu_SO3), [20,20]});

% figure
% obj.PlotShape('b', 0.1,0.1);

%%
scale = 1.2;

sub_obj = EnlargedSuperQuadrics(obj, R, scale);

%% Robot link 6
axes_link = [0.07, 0.18, 0.06];
link = SuperQuadrics({axes_link, [0.50931, 0.40833], [0,0], [0.563725, -0.05271, 0.452293]', [-0.0143401, 0.526693, 0.845612, -0.0856115], [20,20]});
sub_link = EnlargedSuperQuadrics(link, quat2rotm(link.q), scale);

%% Robot link 5
axes_link = [0.13, 0.08, 0.07];
link = SuperQuadrics({axes_link, [0.51867, 1.1162], [0,0], [0.373037, -0.0941539, 0.398832]', [0.188924, 0.0937581, 0.929718, -0.301896], [20,20]});

R_link = zeros(3,3,size(R,3));
for i = 1:size(R,3)
    R_link(:,:,i) = quat2rotm(link.q);
end

sub_link = EnlargedSuperQuadrics(link, R_link, scale);

%% Robot link 4
axes_link4 = [0.15, 0.08, 0.08];
link4 = SuperQuadrics({axes_link4, [0.52535, 1.1094], [0,0], [0.250555, -0.0566034, 0.430801]', [0.91686, -0.363208, -0.114588, 0.119653], [20,20]});
sub_link4 = EnlargedSuperQuadrics(link4, quat2rotm(link4.q), scale);

%% Robot link 3
axes_link3 = [0.15, 0.065, 0.085];
link3 = SuperQuadrics({axes_link3, [0.4983, 0.66919], [0,0], [ 0.0864744, -0.000215062, 0.363642]', [0.0788931, -0.337309, -0.64924, 0.677116], [20,20]});
sub_link3 = EnlargedSuperQuadrics(link3, quat2rotm(link3.q), scale);

%% PCD
figure; axis equal; hold on
link.PlotShape('b',0.3)
obj.PlotShape('g', 0.4)
link4.PlotShape('m',0.3)
link3.PlotShape('b',0.3)

%%
xx = obj.tc - link.tc;
Sigmax = Sigma_t;
 [prob_enlarged, ~, ~, ~, ~] =linearChanceConstraintEnlargedSQ(sub_link, sub_obj, xx, Sigmax, 0)

%%
figure; hold on; axis equal
obj.PlotShape('g', 0.3);
for i=1:size(R,3)
    obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], T(1:3,4,i), rotm2quat(R(:,:,i)), [20,20]});
    obj_i.PlotShape('g', 0.2,0.2);
    pause(0.1);
end
sub_obj.PlotShape('b', 0.6);

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
