clc;close all

%% Loading rotations from samples
T = csv2group(tinychair, 'PCG3');
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

scale = 1.3;

%% tiny chair
axes = [0.2117, 0.187, 0.1331];
tc = [ 0.3993, 0.1287, 0.1911]';
quat = [0.378982, -0.150378, 0.713673, 0.569587];

obj = SuperQuadrics({axes, [1 1], [0,0], tc, quat, [20,20]});
sub_obj = EnlargedSuperQuadrics(obj, R, 1);

%% box1
axes=[0.06,0.3,0.065]';
obj = SuperQuadrics({axes, [0.2,0.2], [0,0], mu_t(1:3,4), rotm2quat(mu_SO3), [20,20]});
sub_obj = EnlargedSuperQuadrics(obj, R, scale);

%% Robot link 1
axes_link = [0.10398,0.10553,0.17304];
link = SuperQuadrics({axes_link, [1 1], [0,0], [-0.036771, 0.00092315, 0.014919]', [0.115937, -0.407734, 0.897055, -0.124916], [20,20]});
sub_link = EnlargedSuperQuadrics(link, link.q, scale);

%% Robot link 9
axes_link = [0.001,0.001,0.001];
link = SuperQuadrics({axes_link, [1 1], [0,0], [0.798997, -0.0442787, 0.358802]', [0.707353, -0.00893546, -0.706743, 0.00927263], [20,20]});
sub_link = EnlargedSuperQuadrics(link, link.q, scale);

%% Robot link 6
axes_link = [0.06723, 0.080807, 0.23949];
link = SuperQuadrics({axes_link, [1 1], [0,0], [0.563762, -0.0523705, 0.452235]', [-0.123494, 0.735943, 0.195778, 0.636245], [20,20]});
sub_link = EnlargedSuperQuadrics(link, quat2rotm(link.q), scale);

%%
xx = obj.tc - link.tc;
% xx = zeros(3,1)
Sigmax = Sigma_t;

%% PCD enlarged
figure; axis equal; hold on
sub_link.PlotShape('b',0.3)
sub_obj.PlotShape('g', 0.4)

prob_enlarged = linearChanceConstraintEnlargedSQ(sub_link, sub_obj, xx, Sigmax, 1)

%% PCD tangent
figure; axis equal; hold on
link.PlotShape('b',0.3)
obj.PlotShape('g', 0.4)

pcd = linearChanceConstraintSQ(link, obj, xx, Sigmax, 'tangent-point-cfc', 1)

%% PCD center
figure; axis equal; hold on
link.PlotShape('b',0.3)
obj.PlotShape('g', 0.4)

figure; axis equal; hold on

plot_ellipse(Sigmax, xx, 'confidence', 1, ...
    'edgecolor', 'm', 'fillcolor', 'm', 'alpha', 0.3)

m1 = link.GetGradientsCanonical();
minkSum = MinkSumClosedForm(link, obj, quat2rotm(link.q), quat2rotm(obj.q));
x_mink = minkSum.GetMinkSumFromGradient(m1);
N = link.N;
X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

mink_sf = surf(X, Y, Z,...
 'EdgeColor', 'm', 'EdgeAlpha', 1,...
 'FaceAlpha', 0.5, 'FaceColor', 'm'); %plot contact space mud yellow


pcd = linearChanceConstraintSQ(link, obj, xx, Sigmax,  'center-point', 1)

%%  compare robot link1
figure; axis equal; hold on
sub_link.PlotShape('g',0.3)
link.PlotShape('b',0.3)



