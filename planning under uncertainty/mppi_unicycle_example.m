function mppi_unicycle_example()
    clc; clear; close all;

    %% Simulation Parameters
    dt          = 0.1;      % Time step
    N           = 20;       % Prediction horizon (number of steps)
    num_samples = 1000;     % Number of sampled trajectories for MPPI
    lambda      = 1.0;      % Temperature parameter for MPPI
    Sigma       = diag([0.2, 0.2]);  % Noise covariance for controls [v; w]
    Q           = diag([1, 1]);      % State (position) error cost
    R           = diag([0.1, 0.1]);  % Control effort cost
    max_steps   = 100;      % Max simulation steps

    % Initial state (x = [px; py; theta])
    x0   = [0; 0; 0];
    x    = x0;

    % Goal position (no strict orientation requirement here)
    goal = [5; 5];

    % Single obstacle (center + radius)
    obs         = [3; 2];
    obs_radius  = 1.0;

    % Nominal (baseline) control sequence for the horizon
    % We can start with zero controls or a naive guess
    u_seq = zeros(2, N);

    % To store the robot trajectory for plotting
    X_traj = x0;

    %% Main Control Loop
    for t = 1:max_steps
        % Use MPPI to find an updated control sequence and the first control
        [u_opt, best_cost, best_traj, u_seq_updated] = ...
            MPPI_step(x, u_seq, dt, N, num_samples, lambda, Sigma, Q, R, goal, obs, obs_radius);

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
    viscircles(obs', obs_radius, 'Color','r','LineWidth',1.5); % obstacle
    plot(goal(1), goal(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','Goal');
    title('MPPI Trajectory for a 2D Unicycle');
    xlabel('x'); ylabel('y');
    legend('Location','best');
end

%% --- MPPI Update Function ---
function [u_opt, best_cost, best_traj, u_seq_updated] = ...
    MPPI_step(x, u_seq, dt, N, num_samples, lambda, Sigma, Q, R, goal, obs, obs_radius)
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
            if norm(x_roll(1:2) - obs) < obs_radius
                stage_cost = stage_cost + 1e10; % large collision penalty
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
