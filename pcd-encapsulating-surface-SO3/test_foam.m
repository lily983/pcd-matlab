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

%% form board
axes=[0.205,0.28,0.0225]';

obj = SuperQuadrics({axes, [0.2,0.2], [0,0], mu_t(1:3,4), rotm2quat(mu_SO3), [20,20]});

figure; axis equal
obj.PlotShape('g', 0.1,0.1);

scale = 1.3;
sub_obj = EnlargedSuperQuadrics(obj, R, scale);

%% Robot link 1
axes_link_1 = [0.10398,0.10553,0.17304];
link_1 = SuperQuadrics({axes_link_1, [1 1], [0,0], [-0.036771, 0.00092315, 0.014919]', [0.707353, -0.00893546, -0.706743, 0.00927263], [20,20]});
sub_link_1 = EnlargedSuperQuadrics(link_1, quat2rotm([0.707353, -0.00893546, -0.706743, 0.00927263]), 1);

%% PCD enlarged
figure; axis equal; hold on
sub_link_1.PlotShape('b',0.3)
sub_obj.PlotShape('g', 0.4)

prob_enlarged = linearChanceConstraintEnlargedSQ(sub_link_1, sub_obj, zeros(3,1), Sigma_t, 1)

%% PCD
figure; axis equal; hold on
link_1.PlotShape('b',0.3)
obj.PlotShape('g', 0.4)

pcd = linearChanceConstraintSQ(link_1, obj, zeros(3,1), Sigma_t, 'tangent-point-cfc', 1)

%%  compare robot link1
figure; axis equal; hold on
sub_link_1.PlotShape('b',0.3)
link_1.PlotShape('g',0.3)



