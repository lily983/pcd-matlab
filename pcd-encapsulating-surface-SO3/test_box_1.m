clc;close all

%% Loading rotations from samples
T = csv2group(marker40new, 'PCG3');
[T_filter, filter_indices] = filterData(T, eye(3));

get_mean_cov(T_filter, 'PCG', 1);

%  print('s1s2positionerror', '-dpng', '-r600')


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

%%
for i=1:size(R,3)
    t(i,:) = T_filter(1:3,4,i);
    quat(i,:) = rotm2quat(R(:,:,i));
    pose_csv(i,:) = [t(i,:), quat(i,2:4), quat(i,1)];
end


%% box
axes=[0.06,0.3,0.065]';

obj = SuperQuadrics({axes, [0.2,0.2], [0,0], mu_t(1:3,4), rotm2quat(mu_SO3), [20,20]});

figure; axis equal; hold on
obj.PlotShape('g', 0.1,0.1);

scale = 1.2;
sub_obj = EnlargedSuperQuadrics(obj, R, scale);
sub_obj.PlotShape('b',0.3)

%% Robot link 1
axes_link = [0.10398,0.10553,0.17304];
link = SuperQuadrics({axes_link, [1 1], [0,0], [-0.036771, 0.00092315, 0.014919]', [0.115937, -0.407734, 0.897055, -0.124916], [20,20]});
sub_link = EnlargedSuperQuadrics(link, link.q, 1);

%% Robot link 9
axes_link = [0.001,0.001,0.001];
link = SuperQuadrics({axes_link, [1 1], [0,0], [0.798997, -0.0442787, 0.358802]', [0.707353, -0.00893546, -0.706743, 0.00927263], [20,20]});
sub_link = EnlargedSuperQuadrics(link, quat2rotm([0.707353, -0.00893546, -0.706743, 0.00927263]), 1);

%% Robot link 6
axes_link = [0.06723, 0.080807, 0.23949];
link = SuperQuadrics({axes_link, [1 1], [0,0], [0.483134, -0.143459, 0.616491]', [0.253533, 0.71257, 0.224211, 0.614569], [20,20]});
sub_link = EnlargedSuperQuadrics(link, quat2rotm(link.q), 1);

%% Robot link 8
axes_link = [0.08, 0.035, 0.05];
link = SuperQuadrics({axes_link, [0.45555, 0.74547], [0,0], [0.493952, -0.328257, 0.584519]', [0.605644, -0.179448, 0.42023, -0.651461], [20,20]});
R_link = zeros(3,3,size(R,3));
for i = 1:size(R,3)
    R_link(:,:,i) = quat2rotm(link.q);
end
sub_link = EnlargedSuperQuadrics(link, R_link, 1);

%% Robot link 7
axes_link = [0.1, 0.055, 0.08];
link = SuperQuadrics({axes_link, [ 0.47578, 0.8109], [0,0], [0.680431, 0.0471307, 0.640214]', [0.0872668, 0.867257, 0.246305, 0.423773], [20,20]});
sub_link = EnlargedSuperQuadrics(link, quat2rotm(link.q), scale);


%%
xx = obj.tc - link.tc;
Sigmax = Sigma_t;

%% PCD enlarged
% figure; axis equal; hold on
% sub_link.PlotShape('b',0.3)
% sub_obj.PlotShape('g', 0.4)

 [prob_enlarged, a, x_mink, ~, ~] =linearChanceConstraintEnlargedSQ(sub_link, sub_obj, xx, Sigmax, 0)
 
 %% Minksum
 psi=[0.927657, -2.05708];
 minkSum = EnlargedMinkSumClosedForm(sub_link, sub_obj);
 m = minkSum.sub1.sq.GetGradientsFromSpherical(psi); n = m ./ norm(m); n_minkSum = minkSum.sub1.R(:,:,1) * n
x_link =  minkSum.sub1.GetPointsFromNormal( n_minkSum )
x_obs =  minkSum.sub1.GetPointsFromNormal( -n_minkSum )
x_mink = minkSum.GetMinkSumFromNormal(n_minkSum)
%% PCD
figure; axis equal; hold on
link.PlotShape('b',0.3)
obj.PlotShape('g', 0.4)

pcd = linearChanceConstraintSQ(link, obj, xx, Sigmax, 'tangent-point-cfc', 0)

 [~, ~, pt_cls, ~] = collision_cfc(link, obj, 'least-squares');
 scatter3(pt_cls.s1(1), pt_cls.s1(2), pt_cls.s1(3), 'r')
 scatter3(pt_cls.s2(1), pt_cls.s2(2), pt_cls.s2(3), 'r')
 
%%  compare robot link1
figure; axis equal; hold on
sub_link.PlotShape('b',0.3)
link.PlotShape('g',0.3)



