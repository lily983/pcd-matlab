clc;close all

%% Loading rotations from samples
T = csv2group(summaryposedata, 'PCG3');
T_filter = filterData(T, eye(3));

R = T_filter(1:3,1:3,:);
mu = get_mean_cov(R, 'SO', 0);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 0);

[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R, 'SO', 0);

%% Loading obj parameters
% axes= rand(3,1) + 0.01*ones(3,1);
% axes=[0.233461976125853;0.322386342019254;0.594523478528556];
% eps=[1.6698    0.5909];
% axes=[0.609438249091412;0.810522765803851;0.115068771291171];
% eps=[1.4116    1.4949];

axes=[0.233461976125853;0.322386342019254;0.594523478528556];
eps=[1.9    0.5909];
sq = SuperQuadrics({axes, eps, [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

figure
sq.PlotShape('b', 0.1,0.1);
%% Scaling constant k
k=1.2;
%% File name
name = "sq-k_1-4";
filetype= "png";

%%
x_new_method_5 = method5(sq, R, k);
visualize_enlarged_surface(sq, R, x_new_method_5)
% filename = name + "_method_5." + filetype;  % change to .pdf if you prefer
% saveas(gcf, filename);
%%
x_new_method_6 = method6(sq, k);
visualize_enlarged_surface(sq, R, x_new_method_6)
% filename = name + "_method_6." + filetype;  % change to .pdf if you prefer
% saveas(gcf, filename);
%%
x_new_method_9 = method9(sq, R, k);
visualize_enlarged_surface(sq, R, x_new_method_9)
% filename = name + "_method_9." + filetype;  % change to .pdf if you prefer
% saveas(gcf, filename);
%%
x_new_method_10 = method10(sq, R, k);
visualize_enlarged_surface(sq, R, x_new_method_10)
% filename = name + "_method_10." + filetype;  % change to .pdf if you prefer
% saveas(gcf, filename);
%%
sq.N=[20,20];
test_sq_10 = EnlargedSuperQuadrics(sq, R, k)
figure
test_sq_10.PlotShape('b', 0.4,0.8)