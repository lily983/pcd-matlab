close all; clear; clc;
rosshutdown;

filePath = '../../data/';
fileName = 'obstacle_';

setenv('ROS_MASTER_URI', 'http://localhost:11311');
rosinit;

%% Setup subscriber and publisher
% Subscribe to clusters of obstacle point cloud
num_pc_sub = rossubscriber('/pc_cluster_number');

shape_pub = rospublisher('/obstacles_config', 'std_msgs/Float32MultiArray');
shape_msg = rosmessage(shape_pub);

pause(1);

while(1)
    disp('Start superellipse generating process...')
    
    num_pc_msg = receive(num_pc_sub);
    num_pc = num_pc_msg.Data;
    
    %% Read from ros topic
    pts = cell(1, num_pc);
    disp(['Totally ', num2str(num_pc), ' clusters'])
    
    try
        for i = 1:num_pc
            disp(['Cluster #: ', num2str(i)])
            
            pause(.1)
            %         pc_sub = rossubscriber(['/obstacles_point_cloud_', num2str(i-1)]);
            %         pause(.2);
            %
            %         pc_msg = receive(pc_sub);
            %
            %         % Read coordinates from PointCloud2
            %         pts{i} = double(readXYZ(pc_msg));
            
            % Read .ply file
            ptsObj = pcread([filePath, fileName, num2str(i-1), '.ply']);
            pts{i} = double(ptsObj.Location);
        end
    catch
        warning('Cannot read PLY file, will try again...')
        continue;
    end
    disp('Retrieved clusters! Start fitting...')
    
    %% Fit superellise
    shapes = nan(num_pc, 6);
    
    for i = 1:num_pc
        pts_fit = [pts{i}(:,1), pts{i}(:,2)];
        
        SQ_par = sq_fit_2d(pts_fit);
        [a0, tc0, th0] = SQ_par.estimate_size;
        x0 = [a0', 1, tc0', th0];
        SQ_par.optimize(x0);
        
        % Save data
        shapes(i,:) = [SQ_par.a(1), SQ_par.a(2), SQ_par.eps,...
            SQ_par.pose.tc(1), SQ_par.pose.tc(2), SQ_par.pose.th];
    end
    
    %% Publish shape data
    shape_msg.Data = shapes';
    send(shape_pub, shape_msg);
    
    csvwrite('../../config/obstacle_config_2D.csv', shapes);
    
    disp('Finished fitting superellipse to each cluster!')
    
end