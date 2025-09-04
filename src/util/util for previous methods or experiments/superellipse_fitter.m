close all; clear; clc;

filePath = '../../data/';
fileName = 'obstacle_';

setenv('ROS_MASTER_URI', 'http://localhost:11311');
rosinit;

% Subscribe to clusters of obstacle point cloud
num_pc_sub = rossubscriber('/pc_cluster_number');
num_pc_msg = receive(num_pc_sub);

%% Get number of clusters
a = dir([filePath, fileName, '*.ply']);
numObj = num_pc_msg.Data;

%% Read and display .ply files
figure; axis equal; grid on; hold on;

for i = 1:numObj
    ptsObj = pcread([filePath, fileName, num2str(i-1), '.ply']);
    pts{i} = double(ptsObj.Location);
    
    plot3(pts{i}(:,1), pts{i}(:,2), pts{i}(:,3),'.')
end

%% Fit an superellipse
figure; axis equal; grid on; hold on;

for i = 1:numObj
    pts_fit{i} = [pts{i}(:,1), pts{i}(:,2)];
    
    SQ_par{i} = sq_fit_2d(pts_fit{i});
    [a0, tc0, th0] = SQ_par{i}.estimate_size;
    x0{i} = [a0', 1, tc0', th0];
    SQ_par{i}.optimize(x0{i});
    
    %% Plots
    SQ_est{i} = SuperEllipse([a0(1), a0(2), 1, tc0(1), tc0(2),...
        th0, 50], 'c', 0);
    
    SQ_fitted{i} = SuperEllipse([SQ_par{i}.a(1), SQ_par{i}.a(2),...
        SQ_par{i}.eps, SQ_par{i}.pose.tc(1), SQ_par{i}.pose.tc(2),...
        SQ_par{i}.pose.th, 50], 'y', 0);
    
    SQ_fitted{i}.PlotShape;
    SQ_est{i}.PlotShape;
    plot(pts_fit{i}(:,1), pts_fit{i}(:,2), '.')
    plot(SQ_par{i}.pt_conv(1,:), SQ_par{i}.pt_conv(2,:), 'm.')
    plot(tc0(1,:), tc0(2,:), 'or')
    
    %% Save data
    shapes(i,:) = [SQ_par{i}.a(1), SQ_par{i}.a(2), SQ_par{i}.eps,...
        SQ_par{i}.pose.tc(1), SQ_par{i}.pose.tc(2), SQ_par{i}.pose.th];
end

csvwrite('../../config/obstacle_config_2D.csv', shapes);

rosshutdown;
