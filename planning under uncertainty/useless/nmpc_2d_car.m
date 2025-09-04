
% NMPC for a 2D Car Kinematic Model with Receding Horizon
% State: [x; y; psi]
% Control: [v; omega]
%
% We re-solve at each time step, applying the first control and
% then updating the state for the next iteration.

clear; clc; close all;

%% Simulation Parameters
dt = 0.05;             % Time step
N = 20;               % Prediction horizon (number of steps in each NMPC solve)
simSteps = 80;        % Total steps in the simulation
x0 = [0; 0; 0];       % Initial state: (x=0, y=0, psi=0)

% Goal position (we only penalize final position, not orientation)
goal = [10; 10];  

% Car safety radius for collision checks
r_car = 0.3;

% Define obstacles (circle centers & radii)
obstacles(1).center = [5; 4];
obstacles(1).radius = 1.0;
obstacles(2).center = [7; 7];
obstacles(2).radius = 1.0;

%% Storage for the trajectory
x_history = zeros(3, simSteps+1);
x_history(:,1) = x0;

% Options for the solver
options = optimoptions('fmincon', 'Display','none', ...
    'Algorithm','sqp', 'MaxFunctionEvaluations',1e5, ...
    'ConstraintTolerance',1e-6, 'OptimalityTolerance',1e-8);

% Initial guess for controls over horizon (v=0, omega=0)
U_init = zeros(2*N,1);

% Receding Horizon Loop
x_current = x0;
for t = 1:simSteps

    % Define cost & constraints for the local horizon
    costFun = @(U) costFunction(U, x_current, goal, N, dt);
    nonlcon = @(U) nonlconFunction(U, x_current, obstacles, r_car, N, dt);

    % Solve the NMPC problem for horizon N
    [U_opt, ~, exitflag, ~] = fmincon(costFun, U_init, [], [], [], [], [], [], nonlcon, options);

    % Check if solver succeeded
    if exitflag <= 0
        warning('Solver did not converge at step %d. Applying zero control.', t);
        u_applied = [0;0];
    else
        % Extract the first control input
        u_applied = U_opt(1:2);
    end

    % Apply the first control to the car
    x_next = carKinematics(x_current, u_applied, dt);
    x_history(:,t+1) = x_next;
    x_current = x_next;

    % Optional: shift initial guess for next time
    % so next iteration starts near last solution
    U_init = [U_opt(3:end); 0; 0];
end

%% Plot Results
figure; hold on; grid on; axis equal;
% Plot trajectory
plot(x_history(1,:), x_history(2,:), 'b-o','LineWidth',2);
% Plot goal
plot(goal(1), goal(2), 'g*','MarkerSize',10);

% Plot obstacles (inflated by r_car)
th = linspace(0,2*pi,100);
for i = 1:length(obstacles)
    ox = obstacles(i).center(1);
    oy = obstacles(i).center(2);
    r_obs = obstacles(i).radius + r_car; % inflate obstacle radius
    plot(ox + r_obs*cos(th), oy + r_obs*sin(th), 'r-','LineWidth',2);
end
xlabel('x (m)'); ylabel('y (m)');
title('2D Car NMPC with Receding Horizon Collision Avoidance');


%% --- Sub-Functions ---

function x_next = carKinematics(x, u, dt)
% Kinematic update for the 2D car
% x = [x; y; psi], u = [v; omega]

    v = u(1);
    omega = u(2);

    x_next = zeros(3,1);
    x_next(1) = x(1) + dt * v * cos(x(3)); % x_{k+1}
    x_next(2) = x(2) + dt * v * sin(x(3)); % y_{k+1}
    x_next(3) = x(3) + dt * omega;        % psi_{k+1}
end

function [J] = costFunction(U, x0, goal, N, dt)
% Cost function for the NMPC problem
% Minimizes final distance to the goal plus small penalty on control usage

    % Simulate the car over the horizon
    [X_traj] = simulateTrajectory(x0, U, N, dt);

    % Terminal position cost: squared distance to goal
    terminal_pos = X_traj(1:2,end);
    dist_goal = norm(terminal_pos - goal)^2;

    % Control effort cost
    U_mat = reshape(U, 2, []);
    ctrl_effort = sum(sum(U_mat.^2));

    % Weights
    w_dist = 100;    % emphasis on reaching goal
    w_ctrl = 0.1;    % emphasis on minimizing control

    % Combine
    J = w_dist * dist_goal + w_ctrl * ctrl_effort;
end

function [c, ceq] = nonlconFunction(U, x0, obstacles, r_car, N, dt)
% Nonlinear constraints for collision avoidance
% c(k) <= 0 for each time step and each obstacle
%    (r_obs + r_car)^2 - ((x - ox)^2 + (y - oy)^2) <= 0

    [X_traj] = simulateTrajectory(x0, U, N, dt);
    numSteps = size(X_traj,2);

    c = [];
    for i = 1:length(obstacles)
        ox = obstacles(i).center(1);
        oy = obstacles(i).center(2);
        r_obs = obstacles(i).radius;

        for k = 1:numSteps
            xk = X_traj(1,k);
            yk = X_traj(2,k);
            % Inequality: (r_obs + r_car)^2 - ( (xk - ox)^2 + (yk - oy)^2 ) <= 0
            c(end+1,1) = (r_obs + r_car)^2 - ((xk - ox)^2 + (yk - oy)^2);
        end
    end
    ceq = [];
end

function [X_traj] = simulateTrajectory(x0, U, N, dt)
% Forward-simulate the 2D car for N steps from initial state x0
% using the control sequence U (length 2*N)

    X_traj = zeros(3, N+1);
    X_traj(:,1) = x0;
    U_mat = reshape(U, 2, []); % 2xN

    for k = 1:N
        x_current = X_traj(:,k);
        u_k = U_mat(:,k);
        X_traj(:,k+1) = carKinematics(x_current, u_k, dt);
    end
end
