% Test script for fitting superquadrics to robot link meshes, using
% probabilistic EMS algorithm
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    convhull, sq_fitting, SuperQuadrics

close all; clear; clc;

addpath('../Superquadrics_Fitting/')
addpath('../Superquadrics_Fitting/functions')

disp('*****************************************************************')
disp('*** Testing for superquadric fitting, using Probabilistic EMS ***')
disp('*****************************************************************')

stl_prefix = 'mesh/';
robot_name = 'panda';
files = dir( [stl_prefix, robot_name, '/*.stl'] );
num_files = size(files, 1);

sq_out_prefix = 'result/';

%% Fitting superquadrics to convex polyhedra
disp('Fitting superquadrics to .stl mesh...')

sq_shape_list = nan(num_files, 12);

vol_poly = nan(1, num_files);
vol_sq = nan(1, num_files);
rel_vol_sq_poly = nan(1, num_files);
error_rel_sq_poly = nan(1, num_files);
error_pt_surf = nan(1, num_files);

for i = 1:num_files
    stl_file_name = [files(i).folder, '/', files(i).name];
    
    % Generate random convex polyhedra
    point_cloud = stlread(stl_file_name).Points';
    [k, vol_poly(i)] = convhull(point_cloud(1,:), point_cloud(2,:),...
        point_cloud(3,:));
    poly = point_cloud(:,unique(k));
    
    % Fit superquadric model, using probabilistic EMS algorithm
    % fitting parameters
    para.w = 0.9;% 0.8 0.7
    para.iterEM = 20; % 20
    para.toleranceEM = 1e-3; % 1e-3
    para.relative_toleranceEM = 0.1; % 0.1
    para.iterLSQ_max = 2; % 2
    para.sigma2 = 0.2; % 0.25
    para.max_switch = 2; % 2
    para.adaptive_upper = 1;
    para.rescale = 0;
    
    para.debug_mode = 0;
    para.arclength = 0.2; %0.2
    
    para.maximum_layer = 3;
    para.minDistance = 0.6; %1
    para.minPoints = 100; %30
    
    [x, point_seg, point_outlier] = SuperquadricsGaussian(poly, para);
    sq_shape.eps = x(1:2);
    sq_shape.a = x(3:5);
    sq_shape.pose.q = eul2quat(x(6:8));
    sq_shape.pose.tc = x(9:end)';
    
    % Save superquadric shape into .csv file
    sq_shape_list(i,:) = [sq_shape.a, sq_shape.eps,...
        sq_shape.pose.tc', quat2axang(sq_shape.pose.q)];
    
    sq = SuperQuadrics({sq_shape.a, sq_shape.eps, [0,0],...
        sq_shape.pose.tc, sq_shape.pose.q, [50, 50]});
    vol_sq(i) = sq.GetVolume();
    
    poly_canonical = quat2rotm(sq.q)\(poly - sq.tc);
    error_pt_surf(i) = sum( abs(sq.GetImplicitFunction(poly_canonical)) )...
        / size(poly_canonical, 2);
    
    % Relative volume
    error_rel_sq_poly(i) = abs(vol_sq(i)-vol_poly(i))./vol_poly(i);
    rel_vol_sq_poly(i) = vol_sq(i)./vol_poly(i);
    
    % Plot link
    figure; hold on; axis equal; axis off;
    lightangle(gca,45,30);
    lighting gouraud;
    
    plot3(point_cloud(1,:), point_cloud(2,:), point_cloud(3,:), 'bo');
    
    plot3(poly(1,:), poly(2,:), poly(3,:), 'r*');
    trisurf(k, point_cloud(1,:), point_cloud(2,:), point_cloud(3,:),...
        'FaceColor', 'y', 'FaceAlpha', 0.7)
    
    sq.PlotShape('g', 0.3);
end

disp(['Mean error avg point to fitted surface distance: ',...
    num2str( mean(error_pt_surf) )]);
disp(['Mean error rel volume SQ fit Polyhedra: ',...
    num2str( mean(error_rel_sq_poly) )]);
disp(['Mean relative volume SQ fit Polyhedra: ',...
    num2str( mean(rel_vol_sq_poly) )]);

csvwrite([sq_out_prefix, robot_name, '_ems.csv'], sq_shape_list);