%% Test Script for Probabilistic Collision Avoidance
% This script tests the probabilistic collision avoidance function
% with different scenarios.

clear all;
close all;
clc;

%% Initialize YALMIP if not already done
if ~exist('sdpvar', 'file')
    error('YALMIP not found. Please install YALMIP to run this code.');
end

%% Common parameters for all tests
params = struct();
params.dt = 0.05;               % Time step (s)
params.N_NMPC = 20;                 % Prediction horizon
params.delta = 0.1;         % Collision probability threshold (10%)
params.R = diag([0.5, 0.1]);     % Stage control cost matrix
% params.Q_k = diag([5, 5]);     % Stage cost matrix
params.Q = diag([50,50]);     % Terminal cost matrix
params.v_min = -0.2;           % Minimum velocity (m/s)
params.v_max = 0.5;            % Maximum velocity (m/s)
params.omega_min = -0.5;       % Minimum angular velocity (rad/s)
params.omega_max = 0.8;        % Maximum angular velocity (rad/s)
params.radius = 0.3;           % Robot radius (m)
params.n_x = 3;
params.n_u = 2;
params.x0 = [0.0; 0.0; 0];  
params.goal = [2; 2; 0];         
direction = params.goal(1:2)-params.x0(1:2);
params.x0(3) = atan(direction(2)/direction(1)); 
params.goal(3) = params.x0(3);
params.goal_tolerance = 0.1;        
params.orientation_tolerance = 0.3;
params.t_final = 2;
params.sim_steps = ceil(params.t_final/params.dt);
params.Sigma_0 = [0.001, 0; 0, 0.001]; % Initial position uncertainty
params.Sigma_propogagtion = diag([0.001, 0.001]);

robot_sq = SuperEllipse([params.radius, params.radius, 1, 0, params.x0(1), params.x0(2), 0.0, 20]);
params.sq = robot_sq;

%%
% Define obstacles
obstacles = {};
obs_sq = SuperEllipse([0.4, 0.2, 1, 0, 1, 1, 0.4, 20]);
obstacles{1} = struct('pos', obs_sq.tc, 'cov', diag([0.001, 0.001]), ...
                     'shape', "ellipse", 'semi_axes', obs_sq.a, 'orientation', obs_sq.ang,...
                     'sq', obs_sq);
% obs_sq_2 = SuperEllipse([0.2, 0.5, 1, 0, 1.5, 0.1, 0.0, 20]);
% obstacles{2} = struct('pos', obs_sq_2.tc, 'cov', diag([0.001, 0.001]), ...
%                      'shape', "ellipse", 'semi_axes', obs_sq_2.a, 'orientation', obs_sq_2.ang,...
%                      'sq', obs_sq_2);
                 
% e_s_plot_trajectory_with_uncertainty(d_tangent_record, obstacles)  
                

