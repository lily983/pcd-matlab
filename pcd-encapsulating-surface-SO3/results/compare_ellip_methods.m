clc; close all
%% Loading rotations from samples
T = csv2group(summaryposedata, 'PCG3');
T_filter = filterData(T, eye(3));

R = T_filter(1:3,1:3,:);
mu = get_mean_cov(R, 'SO', 0);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 0);

[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R, 'SO', 0);

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
ellip = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});

% figure
% ellip.PlotShape('b', 0.1,0.1);
%% Scaling constant k
k=1.4;
%% File name
name = "k_1-4";
filetype= "png";
%% Get encapsulating surface points
x_new_method_1 = method1(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_1)
filename = name + "_method_1." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_2 = method2(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_2)
filename = name + "_method_2." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_3 = method3(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_3)
filename = name + "_method_3." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_4 = method4(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_4)
filename = name + "_method_4." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_5 = method5(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_5)
filename = name + "_method_5." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_6 = method6(ellip, k);
visualize_enlarged_surface(ellip, R, x_new_method_6)
filename = name + "_method_6." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_7 = method7(ellip, R);
visualize_enlarged_surface(ellip, R, x_new_method_7)
filename = name + "_method_7." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_8 = method8(ellip, R);
visualize_enlarged_surface(ellip, R, x_new_method_8)
filename = name + "_method_8." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_9 = method9(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_9)
filename = name + "_method_9." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
%%
x_new_method_10 = method10(ellip, R, k);
visualize_enlarged_surface(ellip, R, x_new_method_10)
filename = name + "_method_10." + filetype;  % change to .pdf if you prefer
saveas(gcf, filename);
