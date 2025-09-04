%% Optimized parameters for NMPC car control
clear all;
close all;
clc;

%% Initialize YALMIP if not already done
if ~exist('sdpvar', 'file')
    error('YALMIP not found. Please install YALMIP to run this code.');
end

%% Common parameters for all tests
params = struct();

% Time and simulation parameters
params.dt = 0.05;               % Time step (s)
params.N_NMPC_base = 8;         % Base prediction horizon (reduced for better performance)
params.N_NMPC_max = 20;         % Maximum prediction horizon 
params.N_NMPC = params.N_NMPC_base; % For compatibility
params.N_NMPC_current = params.N_NMPC_base; % Current horizon (for adaptive horizon)
params.adaptive_horizon = true; % Enable adaptive horizon length
params.t_final = 50;            % Simulation time
params.sim_steps = ceil(params.t_final/params.dt);

% Collision avoidance parameters
params.delta = 0.05;         % Collision probability threshold (5%)
params.radius = 0.3;         % Robot radius (m)

% Cost matrices - adjusted for better balance
params.R = diag([1.0, 2.0]);     % Increased control costs (higher penalty on angular velocity)
params.Q_k = diag([4, 4, 0.5]);  % Stage cost (reduced orientation weight)
params.Q = diag([40, 40, 5]);    % Terminal cost (reduced orientation weight)

% Control limits
params.v_min = 0;              % Minimum velocity (m/s)
params.v_max = 0.3;            % Maximum velocity (m/s)
params.omega_min = -0.2;       % Minimum angular velocity (rad/s)
params.omega_max = 0.2;        % Maximum angular velocity (rad/s)

% Goal parameters
params.x0 = [0.0; 0.0; 0];     % Initial state 
params.goal = [2; 2; 0];       % Goal state

% Calculate initial orientation towards goal
direction = params.goal(1:2)-params.x0(1:2);
params.x0(3) = atan2(direction(2), direction(1)); % Using atan2 for proper quadrant
params.goal(3) = params.x0(3);  % Match goal orientation to initial

params.n_x = 3;
params.n_u = 2;

% Goal reaching parameters
params.goal_tolerance = 0.1;        % Position tolerance
params.orientation_tolerance = 0.3; % Orientation tolerance

% Uncertainty model
params.Sigma_0 = [0.001, 0; 0, 0.001]; % Initial position uncertainty
params.Sigma_propogagtion = diag([0.0005, 0.0005]); % Reduced uncertainty propagation

% Obstacle avoidance parameters - new
params.obstacle_repulsion = 0.5;  % Strength of the repulsion field
params.repulsion_radius = 0.5;    % Radius of influence for repulsion
params.safe_distance = 0.001;       % Distance considered "safe" from obstacles

% Robot model
robot_sq = SuperEllipse([params.radius, params.radius, 1, 0, params.x0(1), params.x0(2), 0.0, 20]);
params.sq = robot_sq;

%% Define obstacles
obstacles = {};

% First obstacle
obs_sq = SuperEllipse([0.4, 0.2, 1, 0, 0.2, 1.2, 0.4, 20]);
obstacles{1} = struct('pos', obs_sq.tc, 'cov', diag([0.001, 0.001]), ...
                     'shape', "ellipse", 'semi_axes', obs_sq.a, 'orientation', obs_sq.ang,...
                     'sq', obs_sq);

% Second obstacle
obs_sq_2 = SuperEllipse([0.2, 0.5, 1, 0, 1.5, 0.1, 0.0, 20]);
obstacles{2} = struct('pos', obs_sq_2.tc, 'cov', diag([0.001, 0.001]), ...
                     'shape', "ellipse", 'semi_axes', obs_sq_2.a, 'orientation', obs_sq_2.ang,...
                     'sq', obs_sq_2);

% Run the optimized NMPC function
result = optimized_CCMPC_2(obstacles, params);

% Visualize results
visualize_results(result, obstacles);