function nmpc_open_world()
% NMPC for a 2D Car in an Open World with Real-Time Sensing of Obstacles
% The car uses a kinematic model:
%   x_{k+1} = x_k + dt*v_k*cos(psi_k)
%   y_{k+1} = y_k + dt*v_k*sin(psi_k)
%   psi_{k+1} = psi_k + dt*omega_k
%
% The car does not know all obstacles in advance but "senses" obstacles
% within a sensing radius. It then uses these in its NMPC optimization.

    clear; clc; close all;

    %% Simulation Parameters
    dt = 0.1;              % time step (seconds)
    N = 10;                % prediction horizon (number of steps per NMPC solve)
    simSteps = 100;        % total simulation iterations (receding horizon steps)
    x0 = [0; 0; 0];        % initial state: [x; y; psi]
    goal = [15; 15];       % goal position
    r_car = 0.3;           % car safety radius
    sensingRadius = 5;     % sensing radius: only obstacles within 5 m are used

    %% Generate Global Obstacles (the car doesn't know all of these a priori)
    numObstacles = 20;
    globalObstacles = struct('center', {}, 'radius', {});
    for i = 1:numObstacles
       % Place obstacles randomly in the region [-5,20] x [-5,20]
       center = [-5 + 25*rand; -5 + 25*rand];
       radius = 0.5 + 0.5*rand;  % radius between 0.5 and 1.0
       globalObstacles(i).center = center;
       globalObstacles(i).radius = radius;
    end

    %% NMPC Setup
    options = optimoptions('fmincon', 'Display','none', ...
        'Algorithm','sqp', 'MaxFunctionEvaluations',1e5, ...
        'ConstraintTolerance',1e-6, 'OptimalityTolerance',1e-8);
    U_init = zeros(2*N,1); % initial guess for control inputs: [v; omega] over N steps

    %% Simulation Loop (Receding Horizon NMPC)
    x_history = zeros(3, simSteps+1);
    x_history(:,1) = x0;
    x_current = x0;
    
    for t = 1:simSteps
       % "Sense" obstacles: extract obstacles within the sensing radius
       localObs = getLocalObstacles(x_current, globalObstacles, sensingRadius);
       
       % Set up NMPC optimization based on the current state and locally sensed obstacles
       costFun = @(U) costFunction(U, x_current, goal, N, dt);
       nonlcon = @(U) nonlconFunction(U, x_current, localObs, r_car, N, dt);
       
       % Solve NMPC problem for horizon N
       [U_opt, ~, exitflag, ~] = fmincon(costFun, U_init, [], [], [], [], [], [], nonlcon, options);
       if exitflag <= 0
           warning('NMPC did not converge at step %d. Applying zero control.', t);
           u_applied = [0; 0];
       else
           u_applied = U_opt(1:2);
       end
       
       % Update the state using the car's kinematics
       x_next = carKinematics(x_current, u_applied, dt);
       x_history(:, t+1) = x_next;
       x_current = x_next;
       
       % Shift initial guess for the next iteration for faster convergence
       U_init = [U_opt(3:end); zeros(2,1)];
    end

    %% Plot the Results
    figure; hold on; grid on; axis equal;
    % Plot the car's trajectory
    plot(x_history(1,:), x_history(2,:), 'b-o', 'LineWidth', 2);
    % Plot the goal
    plot(goal(1), goal(2), 'g*', 'MarkerSize', 10);
    % Plot all global obstacles for reference
    th = linspace(0,2*pi,100);
    for i = 1:length(globalObstacles)
       ox = globalObstacles(i).center(1);
       oy = globalObstacles(i).center(2);
       % Inflate obstacle radius by car safety radius for collision avoidance
       r_obs = globalObstacles(i).radius + r_car;
       plot(ox + r_obs*cos(th), oy + r_obs*sin(th), 'r-', 'LineWidth', 1);
    end
    xlabel('x (m)'); ylabel('y (m)');
    title('NMPC 2D Car in an Open World with Real-Time Sensing');
end

%% --- Sub-Functions ---

function x_next = carKinematics(x, u, dt)
% carKinematics computes the next state of the car.
% x = [x; y; psi], u = [v; omega]
    v = u(1);
    omega = u(2);
    x_next = zeros(3,1);
    x_next(1) = x(1) + dt*v*cos(x(3));
    x_next(2) = x(2) + dt*v*sin(x(3));
    x_next(3) = x(3) + dt*omega;
end

function [X_traj] = simulateTrajectory(x0, U, N, dt)
% simulateTrajectory rolls out the 2D car dynamics for N steps.
    X_traj = zeros(3, N+1);
    X_traj(:,1) = x0;
    U_mat = reshape(U, 2, []);
    for k = 1:N
       X_traj(:, k+1) = carKinematics(X_traj(:,k), U_mat(:,k), dt);
    end
end

function J = costFunction(U, x0, goal, N, dt)
% costFunction computes the NMPC cost.
% It penalizes the squared distance from the final position to the goal
% and a small control effort.
    X_traj = simulateTrajectory(x0, U, N, dt);
    terminal_cost = norm(X_traj(1:2,end) - goal)^2;
    U_mat = reshape(U, 2, []);
    control_cost = sum(sum(U_mat.^2));
    w_terminal = 100;  % weight for terminal cost
    w_control = 0.1;   % weight for control effort
    J = w_terminal*terminal_cost + w_control*control_cost;
end

function [c, ceq] = nonlconFunction(U, x0, obstacles, r_car, N, dt)
% nonlconFunction imposes collision avoidance constraints.
% For each time step and each sensed obstacle, enforce:
%   (r_obs + r_car)^2 - [(x - o_x)^2 + (y - o_y)^2] <= 0.
    X_traj = simulateTrajectory(x0, U, N, dt);
    numSteps = size(X_traj,2);
    c = [];
    for i = 1:length(obstacles)
        ox = obstacles(i).center(1);
        oy = obstacles(i).center(2);
        r_obs = obstacles(i).radius;
        for k = 1:numSteps
            xk = X_traj(1,k);
            yk = X_traj(2,k);
            c(end+1,1) = (r_obs + r_car)^2 - ((xk - ox)^2 + (yk - oy)^2);
        end
    end
    ceq = [];
end

function localObs = getLocalObstacles(x, globalObs, sensingRadius)
% getLocalObstacles returns obstacles from globalObs that are within
% the sensingRadius of the car's current position x.
    localObs = [];
    for i = 1:length(globalObs)
       center = globalObs(i).center;
       if norm(x(1:2) - center) <= sensingRadius
           localObs = [localObs, globalObs(i)];
       end
    end
end
