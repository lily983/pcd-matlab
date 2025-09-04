
    clc; clear; close all;
% 
%     %% Simulation Parameters
    dt          = 0.1;      % Time step
    N           = 20;       % Prediction horizon (number of steps)
    num_samples = 1000;     % Number of sampled trajectories for MPPI
    lambda      = 1.0;      % Temperature parameter for MPPI
    Sigma       = diag([0.1, 0.1]);  % Noise covariance for controls [v; w]
    Q           = diag([1, 1]);      % State (position) error cost
    R           = diag([0.0, 0.0]);  % Control effort cost
    Q_N           = diag([20, 20]);  % Control effort cost
    max_steps   = 200;      % Max simulation steps
    
    Sigma_pos = diag([0.005,0.005]);
    threshold = 0.1;

    % Initial state (x = [px; py; theta])
    x0   = [0; 1; 0];
    x    = x0;

    % Goal position (no strict orientation requirement here)
    goal = [5; -1];

    % Define robot parameters (for the super-ellipse representing the robot)
    robot_sq_param = struct();
    robot_sq_param.a    = [0.6, 0.6];
    robot_sq_param.eps  = 1.0;
    robot_sq_param.tc   = x0(1:2);
    robot_sq_param.ang  = 0.0;
    
    % Define multiple obstacles as an array of structures.
    % Each obstacle has fields: a, eps, tc (center), and ang.
    obs1 = struct('a', [0.5, 0.8], 'eps', 1.0, 'tc', [4; 2], 'ang', 0.2);
    obs2 = struct('a', [0.4, 0.6], 'eps', 1.0, 'tc', [1.8; -0.4], 'ang', -0.1);
    % You can add more obstacles as needed:
    % obs3 = struct('a', [...], 'eps', ..., 'tc', [...], 'ang', ...);
    obs_sq_params = [obs1, obs2];

    % Nominal (baseline) control sequence for the horizon
    % We can start with zero controls or a naive guess
    u_seq = zeros(2, N);

    % To store the robot trajectory for plotting
    X_traj = x0;

    %% Main Control Loop
    tic
    for t = 1:max_steps
        % Use MPPI to find an updated control sequence and the first control.
        % Note that we now pass the entire array of obstacle parameters.
        [u_opt, best_cost, best_traj, u_seq_updated] = ...
            MPPI_step(x, u_seq, dt, N, num_samples, lambda, Sigma, Q, R, goal, obs_sq_params, robot_sq_param, Sigma_pos, threshold, Q_N);

        % Apply the first control from the updated sequence
        x = unicycle_dynamics(x, u_opt, dt);
        X_traj = [X_traj, x];

        % Shift the horizon: discard the first control and append a new one
        u_seq = [u_seq_updated(:,2:end), u_seq_updated(:,end)];

        % Stop if we are close enough to the goal
        if norm(x(1:2) - goal) < 0.2
            disp('Goal reached (within tolerance)!');
            break;
        end
    end
    t_ellip = toc;
    %% Plot Results
    figure; hold on; axis equal;
    plot(X_traj(1,:), X_traj(2,:), 'b-o', 'LineWidth', 1.5, 'DisplayName','Robot Trajectory');
%     visualizeEllipseRobot(X_traj', robot_sq_param.a(1), robot_sq_param.a(2));
    
    % Loop through each obstacle to visualize them
    for i = 1:length(obs_sq_params)
        visualizeEllipse(obs_sq_params(i).a(1), obs_sq_params(i).a(2), obs_sq_params(i).tc', obs_sq_params(i).ang);
        visualize_sq_error(obs_sq_params(i), Sigma_pos)
    end
    plot(goal(1), goal(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','Goal');
%     title('MPPI Trajectory for a 2D Unicycle');
%     xlabel('x'); ylabel('y');
%     legend('Location','best');


%% --- MPPI Update Function ---
function [u_opt, best_cost, best_traj, u_seq_updated] = ...
    MPPI_step(x, u_seq, dt, N, num_samples, lambda, Sigma, Q, R, goal, obs_sq_params, robot_sq_param, Sigma_pos, threshold, Q_N)
% x: current state [px; py; theta]
% u_seq: nominal control sequence, size(2,N) = [v; w] for each of N steps
% dt, N: time step and horizon length
% num_samples: number of sample rollouts
% lambda: temperature parameter
% Sigma: 2x2 noise covariance for [v; w]
% Q, R: cost weights
% goal: desired position (2D)
% obs_sq_params: array of obstacle parameters (each a struct with fields 'a', 'eps', 'tc', 'ang')
% robot_sq_param: robot parameters (struct with fields 'a', 'eps', 'tc', 'ang')
% Sigma_pos: position noise covariance (for chance constraints)
% threshold: chance constraint threshold

    costs          = zeros(num_samples, 1);
    control_noise  = zeros(2, N, num_samples);

    % Cholesky factor for sampling from Sigma
    L = chol(Sigma, 'lower');

    best_cost  = inf;
    best_index = 1;
    
    R_robot = angle2rotm(robot_sq_param.ang);
    
    % Sample and evaluate each rollout
    for k = 1:num_samples
        % Sample noise for the entire horizon
        noise_k = L * randn(2, N);

        % Simulate forward from current state
        x_roll     = x;
        total_cost = 0;

        for i = 1:N
            % Control + noise
            u_noisy = u_seq(:, i) + noise_k(:, i);

            % Propagate dynamics
            x_roll = unicycle_dynamics(x_roll, u_noisy, dt);

            % Compute stage cost from position error and control effort
            pos_err    = x_roll(1:2) - goal;
            stage_cost = pos_err' * Q * pos_err + u_noisy' * R * u_noisy;

            % Initialize obstacle penalty for this step
            obstacle_penalty = 0;
            % Loop over each obstacle
            for j = 1:length(obs_sq_params)
                obs_param = obs_sq_params(j);
                % Compute the vector from the obstacle center to the robot
                p_io = x_roll(1:2) - obs_param.tc;
                if norm(p_io) == 0
                    % Avoid division by zero (if robot exactly at obstacle center)
                    continue;
                end
                R_obs = angle2rotm(obs_param.ang);

                S_io = diag(1./(obs_param.a + robot_sq_param.a).^2);
                Omega_io = R_obs'*S_io*R_obs;
                Omega_sqrt = sqrtm(Omega_io); % Square root of matrix

                % In transformed space, the boundary is at distance 1
                b_io = 1;
                
                a_io = p_io./norm(p_io);
                
                % Compute the chance constraint value
                constraint = 0.5 + 0.5 * erf((b_io - a_io' * Omega_sqrt * p_io) /...
                    sqrt(max(2 * a_io' * Omega_sqrt * ((i+1) * Sigma_pos) * Omega_sqrt' * a_io, 1e-6)));
                
                if threshold < constraint
                    % Add a large penalty if the constraint is violated
                    obstacle_penalty = obstacle_penalty + ( - threshold + constraint + 1e4);
                end
            end

            % Add obstacle penalty (if any) to the stage cost
            stage_cost = stage_cost + obstacle_penalty;

%             stage_cost = stage_cost;

            total_cost = total_cost + stage_cost;
        end
        
        terminal_pos    = x_roll(1:2) - goal;
        terminal_cost = terminal_pos' * Q_N * terminal_pos;
            
        total_cost = total_cost + terminal_cost;

        costs(k)             = total_cost;
        control_noise(:,:,k) = noise_k;

        % Track best rollout for debugging
        if total_cost < best_cost
            best_cost  = total_cost;
            best_index = k;
        end
    end

    % Compute path integral weights
    beta    = min(costs);  % for numerical stability
    weights = exp(-(costs - beta) / lambda);
    Z       = sum(weights);

    % Update the nominal control sequence with weighted noise
    u_seq_updated = zeros(size(u_seq));
    for i = 1:N
        delta_u = [0; 0];
        for k = 1:num_samples
            delta_u = delta_u + (weights(k) * control_noise(:, i, k));
        end
        u_seq_updated(:, i) = u_seq(:, i) + (delta_u / Z);
    end

    % Pick the first control as the optimal to apply now
    u_opt = u_seq_updated(:, 1);

    % For debugging: reconstruct the best trajectory
    x_roll   = x;
    best_traj = x_roll;
    best_noise = control_noise(:,:,best_index);
    for i = 1:N
        u_noisy  = u_seq(:, i) + best_noise(:, i);
        x_roll   = unicycle_dynamics(x_roll, u_noisy, dt);
        best_traj = [best_traj, x_roll];
    end
end

%% --- Unicycle Dynamics ---
function x_next = unicycle_dynamics(x, u, dt)
% x = [px; py; theta]
% u = [v; w] (linear velocity, angular velocity)
% dt: time step

    px   = x(1);
    py   = x(2);
    th   = x(3);
    v    = u(1);
    w    = u(2);

    px_next = px + dt * v * cos(th);
    py_next = py + dt * v * sin(th);
    th_next = th + dt * w;
    

    x_next = [px_next; py_next; th_next];
end

%% --- SuperEllipse Functions ---
% These functions compute points on a superellipse given a direction vector.
function x = get_points_from_normal(sq_param, n)
    m = get_gradients_from_normal(sq_param, n);
    x = get_points_from_gradient(sq_param, m);
end

function m = get_gradients_from_normal(sq_param, n)
    dual_sq.a   = 2./(sq_param.a * sq_param.eps);
    dual_sq.eps = 2 - sq_param.eps;
    dual_implicit = get_implicit_function(dual_sq, n);
    m = n./eps_fun(dual_implicit+1, dual_sq.eps/2);
end

function x = get_points_from_gradient(sq_param, m)
    sq_eps = sq_param.eps;
    x(1,:) = sq_param.a(1) * eps_fun(sq_eps * sq_param.a(1)/2 * m(1,:), sq_eps/(2-sq_eps));
    x(2,:) = sq_param.a(2) * eps_fun(sq_eps * sq_param.a(2)/2 * m(2,:), sq_eps/(2-sq_eps));
end

function f = get_implicit_function(sq_param, x)
    f = eps_fun( (x(1,:) ./ sq_param.a(1)).^2, 1/sq_param.eps ) +...
        eps_fun( (x(2,:) ./ sq_param.a(2)).^2, 1/sq_param.eps ) - 1;
end

function val = eps_fun(x, epsilon)
% eps_fun: generalized exponentiation function
    val = sign(x) .* abs(x).^epsilon;
end

%% --- Visualization Functions ---
function visualizeEllipseRobot(traj, a, b)
% VISUALIZEELLIPSEROBOT Visualizes the 2D car (ellipse robot) along a trajectory.
%
%   traj - A Nx3 matrix [x, y, theta] (position and orientation)
%   a, b - semi-major and semi-minor axes lengths of the ellipse

    % Plot the trajectory as a dashed line for context
    plot(traj(:,1), traj(:,2), 'k--', 'LineWidth', 1.5);
    
    % Define the ellipse in its local coordinate system
    t = linspace(0, 2*pi, 50);
    ellipse_local = [a * cos(t); b * sin(t)];
    
    % Loop through each pose in the trajectory
    for i = 1:size(traj, 1)
        x = traj(i, 1);
        y = traj(i, 2);
        theta = traj(i, 3);
        
        % Rotate and translate the ellipse into global coordinates
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        ellipse_global = R * ellipse_local;
        ellipse_global(1,:) = ellipse_global(1,:) + x;
        ellipse_global(2,:) = ellipse_global(2,:) + y;
        
        % Draw the ellipse with transparency
        patch(ellipse_global(1,:), ellipse_global(2,:), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
        plot(x, y, 'ko', 'MarkerFaceColor', 'k');
        pause(0.05);
    end
end

% Assume visualizeEllipse and visualize_sq_error are defined similarly for obstacles.
% For example, visualizeEllipse might look like:
function visualizeEllipse(a, b, tc, ang)
    % Create the ellipse in the local frame
    t = linspace(0, 2*pi, 50);
    ellipse_local = [a * cos(t); b * sin(t)];
    % Rotation matrix based on obstacle orientation
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    ellipse_global = R * ellipse_local;
    ellipse_global(1,:) = ellipse_global(1,:) + tc(1);
    ellipse_global(2,:) = ellipse_global(2,:) + tc(2);
    patch(ellipse_global(1,:), ellipse_global(2,:), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
end

function visualize_sq_error(sq_param, Sigma_pos)
    % This function can visualize the error bounds (or chance constraint region)
    % for the given super-ellipse obstacle. Here, we simply plot the center.
    plot(sq_param.tc(1), sq_param.tc(2), 'mx', 'MarkerSize', 10, 'LineWidth', 2);
end

% --- Helper function for rotation matrix from angle ---
function R = angle2rotm(theta)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end
