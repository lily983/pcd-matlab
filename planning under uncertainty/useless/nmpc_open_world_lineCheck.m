function nmpc_open_world_lineCheck()
% NMPC for a 2D Car in an Open World with Real-Time Sensing and
% additional edge collision checking (sampling along the connecting line)
%
% Kinematic model:
%   x_{k+1} = x_k + dt*v_k*cos(psi_k)
%   y_{k+1} = y_k + dt*v_k*sin(psi_k)
%   psi_{k+1} = psi_k + dt*omega_k
%
% Collision avoidance: For each obstacle (circle) with center (ox,oy)
% and radius r_obs, ensure that the car (with safety radius r_car) remains
% outside the inflated obstacle:
%
%   (r_obs + r_car)^2 - ((x - ox)^2 + (y - oy)^2) <= 0.
%
% In addition to checking the discrete trajectory, we sample points
% along the line segments between consecutive states.

    clear; clc; close all;

    %% Simulation Parameters
    dt = 0.1;              % Time step (seconds)
    N = 10;                % Prediction horizon (steps per NMPC solve)
    simSteps = 100;        % Total number of receding horizon iterations
    x0 = [0; 0; 0];        % Initial state: [x; y; psi]
    goal = [15; 15];       % Goal position (only position is penalized)
    r_car = 0.3;           % Safety radius of the car
    sensingRadius = 5;     % Sensing radius: only obstacles within this are used

    %% Global Obstacles (unknown a priori)
    numObstacles = 20;
    globalObstacles = struct('center', {}, 'radius', {});
    for i = 1:numObstacles
       % Randomly place obstacles in region [-5,20] x [-5,20]
       center = [-5 + 25*rand; -5 + 25*rand];
       radius = 0.5 + 0.5*rand;  % radius between 0.5 and 1.0
       globalObstacles(i).center = center;
       globalObstacles(i).radius = radius;
    end

    %% NMPC Setup
    options = optimoptions('fmincon', 'Display','none', ...
        'Algorithm','sqp', 'MaxFunctionEvaluations',1e5, ...
        'ConstraintTolerance',1e-6, 'OptimalityTolerance',1e-8);
    U_init = zeros(2*N,1); % initial guess for control inputs over horizon

    %% Simulation Loop (Receding Horizon NMPC)
    x_history = zeros(3, simSteps+1);
    x_history(:,1) = x0;
    x_current = x0;
    
    for t = 1:simSteps
       % "Sense" obstacles: only consider those within the sensing radius
       localObs = getLocalObstacles(x_current, globalObstacles, sensingRadius);
       
       % Define NMPC optimization functions with current state and local obstacles
       costFun = @(U) costFunction(U, x_current, goal, N, dt);
       nonlcon = @(U) nonlconFunction(U, x_current, localObs, r_car, N, dt);
       
       % Solve NMPC for horizon N
       [U_opt, ~, exitflag, ~] = fmincon(costFun, U_init, [], [], [], [], [], [], nonlcon, options);
       if exitflag <= 0
           warning('NMPC did not converge at step %d. Applying zero control.', t);
           u_applied = [0; 0];
       else
           u_applied = U_opt(1:2);
       end
       
       % Update state using the kinematic model
       x_next = carKinematics(x_current, u_applied, dt);
       x_history(:, t+1) = x_next;
       x_current = x_next;
       
       % Shift initial guess for faster convergence next iteration
       U_init = [U_opt(3:end); zeros(2,1)];
    end

    %% Plot Results
    figure; hold on; grid on; axis equal;
    % Plot car trajectory
    plot(x_history(1,:), x_history(2,:), 'b-o','LineWidth', 2);
    % Plot goal
    plot(goal(1), goal(2), 'g*','MarkerSize', 10);
    % Plot all global obstacles (inflated by r_car)
    th = linspace(0,2*pi,100);
    for i = 1:length(globalObstacles)
       ox = globalObstacles(i).center(1);
       oy = globalObstacles(i).center(2);
       r_obs = globalObstacles(i).radius + r_car;
       plot(ox + r_obs*cos(th), oy + r_obs*sin(th), 'r-', 'LineWidth', 1);
    end
    xlabel('x (m)'); ylabel('y (m)');
    title('2D Car NMPC with Real-Time Sensing and Edge Collision Check');
end

%% --- Sub-Functions ---

function x_next = carKinematics(x, u, dt)
% Kinematic update for the 2D car
% x = [x; y; psi], u = [v; omega]
    v = u(1);
    omega = u(2);
    x_next = zeros(3,1);
    x_next(1) = x(1) + dt*v*cos(x(3));
    x_next(2) = x(2) + dt*v*sin(x(3));
    x_next(3) = x(3) + dt*omega;
end

function [X_traj] = simulateTrajectory(x0, U, N, dt)
% Forward-simulate the 2D car for N steps given a control sequence U.
    X_traj = zeros(3, N+1);
    X_traj(:,1) = x0;
    U_mat = reshape(U, 2, []);
    for k = 1:N
       X_traj(:, k+1) = carKinematics(X_traj(:,k), U_mat(:,k), dt);
    end
end

function J = costFunction(U, x0, goal, N, dt)
% Cost function for NMPC: penalize terminal distance to goal and control effort.
    X_traj = simulateTrajectory(x0, U, N, dt);
    terminal_cost = norm(X_traj(1:2,end) - goal)^2;
    U_mat = reshape(U, 2, []);
    control_cost = sum(sum(U_mat.^2));
    w_terminal = 100;  % weight for terminal position cost
    w_control = 0.1;   % weight for control effort
    J = w_terminal*terminal_cost + w_control*control_cost;
end

function [c, ceq] = nonlconFunction(U, x0, obstacles, r_car, N, dt)
% Nonlinear constraints: ensure collision avoidance along the trajectory.
% For each obstacle, for each edge (between consecutive states),
% sample intermediate points and enforce:
%   (r_obs + r_car)^2 - ((x - ox)^2 + (y - oy)^2) <= 0.
    X_traj = simulateTrajectory(x0, U, N, dt);
    numSteps = size(X_traj,2);
    numEdgeSamples = 5;  % number of samples between consecutive states
    c = [];
    for i = 1:length(obstacles)
        ox = obstacles(i).center(1);
        oy = obstacles(i).center(2);
        r_obs = obstacles(i).radius;
        % For each edge in the trajectory
        for k = 1:(numSteps-1)
            % Sample additional points along the line between X(:,k) and X(:,k+1)
            for s = 0:numEdgeSamples
                alpha = s / numEdgeSamples;
                pos = (1 - alpha)*X_traj(1:2,k) + alpha*X_traj(1:2,k+1);
                % Constraint: (r_obs + r_car)^2 - ((pos(1)-ox)^2 + (pos(2)-oy)^2) <= 0
                c(end+1,1) = (r_obs + r_car)^2 - ((pos(1)-ox)^2 + (pos(2)-oy)^2);
            end
        end
    end
    ceq = [];
end

function localObs = getLocalObstacles(x, globalObs, sensingRadius)
% Returns obstacles from globalObs that are within the sensingRadius of x.
    localObs = [];
    for i = 1:length(globalObs)
       center = globalObs(i).center;
       if norm(x(1:2) - center) <= sensingRadius
           localObs = [localObs, globalObs(i)];
       end
    end
end
