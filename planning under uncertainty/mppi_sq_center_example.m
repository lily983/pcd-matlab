function mppi_sq_center_example()
    clc; clear; close all;

    %% Simulation Parameters
    dt          = 0.1;      % Time step
    N           = 30;       % Prediction horizon (number of steps)
    num_samples = 1000;     % Number of sampled trajectories for MPPI
    lambda      = 1.0;      % Temperature parameter for MPPI
    Sigma       = diag([0.2, 0.2]);  % Noise covariance for controls [v; w]
    Q           = diag([1, 1]);      % State (position) error cost
    R           = diag([0.1, 0.1]);  % Control effort cost
    max_steps   = 200;      % Max simulation steps
    
    Sigma_pos = diag([0.005,0.005]);
    threshold = 0.05;

    % Initial state (x = [px; py; theta])
    x0   = [0; 1; 0];
    x    = x0;

    % Goal position (no strict orientation requirement here)
    goal = [5; 1];

    % Single obstacle (center + radius)
    obs         = [3; 2];
    obs_radius  = 1.0;
    
    robot_sq_param = struct();
    robot_sq_param.a = [0.8, 0.8];
    robot_sq_param.eps = 1.0;
    robot_sq_param.tc = x0(1:2);
    robot_sq_param.ang = 0.0;
    
    obs_sq_param = struct();
    obs_sq_param.a = [0.5, 0.8];
    obs_sq_param.eps = 1.0;
    obs_sq_param.tc = obs;
    obs_sq_param.ang = 0.2;

    % Nominal (baseline) control sequence for the horizon
    % We can start with zero controls or a naive guess
    u_seq = zeros(2, N);

    % To store the robot trajectory for plotting
    X_traj = x0;

    %% Main Control Loop
    for t = 1:max_steps
        % Use MPPI to find an updated control sequence and the first control
        [u_opt, best_cost, best_traj, u_seq_updated] = ...
            MPPI_step(x, u_seq, dt, N, num_samples, lambda, Sigma, Q, R, goal, obs, obs_sq_param, robot_sq_param, Sigma_pos, threshold);

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

    %% Plot Results
    figure; hold on; axis equal;
    plot(X_traj(1,:), X_traj(2,:), 'b-o', 'LineWidth', 1.5, 'DisplayName','Robot Trajectory');
    visualizeEllipseRobot(X_traj', robot_sq_param.a(1), robot_sq_param.a(2));
%     viscircles(obs', obs_radius, 'Color','r','LineWidth',1.5); % obstacle
    visualizeEllipse(obs_sq_param.a(1), obs_sq_param.a(2), obs_sq_param.tc', obs_sq_param.ang);
    visualize_sq_error(obs_sq_param, Sigma_pos)
    plot(goal(1), goal(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','Goal');
    title('MPPI Trajectory for a 2D Unicycle');
    xlabel('x'); ylabel('y');
    legend('Location','best');
end

%% --- MPPI Update Function ---
function [u_opt, best_cost, best_traj, u_seq_updated] = ...
    MPPI_step(x, u_seq, dt, N, num_samples, lambda, Sigma, Q, R, goal, obs, obs_sq_param, robot_sq_param, Sigma_pos, threshold)
% x: current state [px; py; theta]
% u_seq: nominal control sequence, size(2,N) = [v; w] for each of N steps
% dt, N: time step and horizon length
% num_samples: number of sample rollouts
% lambda: temperature parameter
% Sigma: 2x2 noise covariance for [v; w]
% Q, R: cost weights
% goal: desired position (2D)
% obs, obs_radius: obstacle center and radius

    costs          = zeros(num_samples, 1);
    control_noise  = zeros(2, N, num_samples);

    % Cholesky factor for sampling from Sigma
    L = chol(Sigma, 'lower');

    best_cost  = inf;
    best_index = 1;
    
    R_robot = angle2rotm(robot_sq_param.ang);
    R_obs = angle2rotm(obs_sq_param.ang);

    % Sample and evaluate each rollout
    for k = 1:num_samples
        % Sample noise for the entire horizon
        noise_k = L * randn(2, N);

        % Simulate forward
        x_roll     = x;
        total_cost = 0;

        for i = 1:N
            % Control + noise
            u_noisy = u_seq(:, i) + noise_k(:, i);

            % Optionally clamp controls here if needed
            % e.g. u_noisy(1) = max(min(u_noisy(1), v_max), v_min);

            % Propagate dynamics
            x_roll = unicycle_dynamics(x_roll, u_noisy, dt);

            % Compute stage cost
            %  - position error cost
            pos_err    = x_roll(1:2) - goal;
            stage_cost = pos_err' * Q * pos_err + u_noisy' * R * u_noisy;

            % Obstacle penalty (simple "infinite" penalty if inside obstacle)
            %%% chance constraint 
            p_io = x_roll(1:2) - obs;
            a_io = p_io./norm(p_io);
           
            x_robot = R_robot * get_points_from_normal(robot_sq_param, R_robot' *  a_io);
            x_obs = R_obs * get_points_from_normal(obs_sq_param, R_obs' *  a_io);
            x_mink = x_robot + x_obs;
            
            b_io = norm(x_mink);
            
            constraint = 0.5 + 0.5*erf((b_io - a_io' * p_io)/sqrt(max( 2*a_io'*((i+1) * Sigma_pos)*a_io, 1e-6)));
            
            if threshold < constraint
                stage_cost  = threshold - constraint + 1e8;
            end

            total_cost = total_cost + stage_cost;
        end

        costs(k)             = total_cost;
        control_noise(:,:,k) = noise_k;

        % Track best rollout for debugging
        if total_cost < best_cost
            best_cost  = total_cost;
            best_index = k;
        end
    end

    % Path integral weights
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

%%  --- SuperEllipse ---
% sq_param = {'a', [a, b], 'epsilon', eps, 'tc', [x;y], 'ang', theta}
function x = get_points_from_normal(sq_param, n)

m = get_gradients_from_normal(sq_param, n);

x = get_points_from_gradient(sq_param, m);

end

function m = get_gradients_from_normal(sq_param, n)

dual_sq = struct();
dual_sq.a = 2./(sq_param.a * sq_param.eps);
dual_sq.eps = 2 - sq_param.eps;
% dual_sq.pos = sq_param.tc;
% dual_sq.ang = sq_param.ang;
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
%eps_fun exponentiation function
val = sign(x).*abs(x).^epsilon;
end

%% --- Visualization ---
function visualizeEllipseRobot(traj, a, b)
% VISUALIZEELLIPSEROBOT Visualizes a 2D car (ellipse robot) along a planned trajectory.
%
%   visualizeEllipseRobot(traj, a, b) plots an ellipse representing the robot at 
%   each pose along the trajectory.
%
%   Input:
%       traj - A Nx3 matrix where each row is [x, y, theta] (position and orientation)
%       a    - Semi-major axis length of the ellipse
%       b    - Semi-minor axis length of the ellipse
%
%   Example:
%       traj = [0 0 0; 1 0.5 pi/6; 2 1 pi/4; 3 1.5 pi/3];
%       visualizeEllipseRobot(traj, 0.5, 0.3);
%

    % Open a new figure and configure the plot
%     figure;
%     hold on;
%     grid on;
%     axis equal;
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     title('Robot Trajectory with Ellipse Visualization');
    
    % Plot the entire trajectory path as a dashed line for context
    plot(traj(:,1), traj(:,2), 'k--', 'LineWidth', 1.5);
    
    % Create ellipse in its local coordinate system
    t = linspace(0, 2*pi, 50);
    ellipse_local = [a * cos(t); b * sin(t)];
    
    % Loop through each pose in the trajectory
    for i = 1:size(traj, 1)
        % Extract the pose: position (x,y) and orientation (theta)
        x = traj(i, 1);
        y = traj(i, 2);
        theta = traj(i, 3);
        
        % Build the rotation matrix for current heading
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        % Rotate and translate the ellipse to global coordinates
        ellipse_global = R * ellipse_local;
        ellipse_global(1,:) = ellipse_global(1,:) + x;
        ellipse_global(2,:) = ellipse_global(2,:) + y;
        
        % Draw the ellipse using patch with some transparency
        h = patch(ellipse_global(1,:), ellipse_global(2,:), 'b', ...
            'FaceAlpha', 0.3, 'EdgeColor', 'r');
        
        % Optionally plot the center of the ellipse
        plot(x, y, 'ko', 'MarkerFaceColor', 'k');
        
        % Pause to animate the movement (adjust the delay as needed)
        pause(0.05);
        
        % Uncomment the following line if you want to clear the previous ellipse for animation
        % if i < size(traj, 1), delete(h); end
    end
%     hold off;
end

