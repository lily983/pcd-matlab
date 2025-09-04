%% Minimal NMPC implementation to get the car moving
% This is stripped down to ensure it works, then you can add features back

function d = minimal_NMPC(obstacles, params)
d.p = params;
d.c.obstacles = obstacles;

% Pre-allocate memory
d.s.x = zeros(d.p.n_x, d.p.sim_steps+1);
d.s.u = zeros(d.p.n_u, d.p.sim_steps);
d.s.CPU_time = zeros(d.p.sim_steps, 1);

% Set initial state
d.s.x(:,1) = d.p.x0;

% Directly control the first step to ensure movement
disp('Setting initial movement...');
d.s.u(:,1) = [0.2; 0]; % Small forward velocity, no turning
d.s.x(:,2) = car_dynamics(d.s.x(:,1), d.s.u(:,1), d.p.dt);
disp('Initial state:');
disp(d.s.x(:,1));
disp('After first step:');
disp(d.s.x(:,2));

% Main simulation loop
for t = 2:d.p.sim_steps
    fprintf('Time step %d, Current position: (%.2f, %.2f)\n', t, d.s.x(1,t), d.s.x(2,t));
    
    % Calculate distance to goal
    distance_to_goal = norm(d.s.x(1:2,t) - d.p.goal(1:2));
    if distance_to_goal <= d.p.goal_tolerance
        fprintf('Goal reached at timestep %d!\n', t);
        break;
    end
    
    % Solve simple NMPC problem
    [d, u_next] = simple_NMPC_step(d, t);
    
    % Apply control and evolve dynamics
    d.s.u(:,t) = u_next;
    d.s.x(:,t+1) = car_dynamics(d.s.x(:,t), d.s.u(:,t), d.p.dt);
    
    % Debug: Verify movement
    movement = norm(d.s.x(1:2,t+1) - d.s.x(1:2,t));
    fprintf('Movement: %.4f meters\n', movement);
    if movement < 1e-6
        warning('No movement detected at time %d!', t);
        
        % Force movement if stuck
        direction = d.p.goal(1:2) - d.s.x(1:2,t);
        angle = atan2(direction(2), direction(1));
        d.s.u(:,t) = [0.2; 0.05*sign(angle - d.s.x(3,t))];
        d.s.x(:,t+1) = car_dynamics(d.s.x(:,t), d.s.u(:,t), d.p.dt);
        fprintf('Forced movement to: (%.2f, %.2f)\n', d.s.x(1,t+1), d.s.x(2,t+1));
    end
end

% Display results
disp('Final state:');
disp(d.s.x(:,end));
end

function [d, u_next] = simple_NMPC_step(d, t)
% Simplified NMPC step - no complex constraints

% Parameters
N = 5; % Short horizon for simplicity
current_state = d.s.x(:,t);
obstacles = d.c.obstacles;

tic;

% Formulate the optimization problem
try
    % Variables: U = [u1, u2, ..., uN]
    n_vars = d.p.n_u * N;
    
    % Initial guess - straight line toward goal
    direction = d.p.goal(1:2) - current_state(1:2);
    angle_to_goal = atan2(direction(2), direction(1));
    angle_diff = angle_to_goal - current_state(3);
    angle_diff = mod(angle_diff + pi, 2*pi) - pi; % Normalize to [-pi, pi]
    
    x0 = zeros(n_vars, 1);
    for i = 1:N
        idx = (i-1)*d.p.n_u + 1;
        x0(idx) = 0.2; % Velocity
        x0(idx+1) = 0.05 * sign(angle_diff); % Small turn toward goal
    end
    
    % Lower and upper bounds
    lb = repmat([d.p.v_min; d.p.omega_min], N, 1);
    ub = repmat([d.p.v_max; d.p.omega_max], N, 1);
    
    % Define options
    options = optimoptions('fmincon', ...
                           'Algorithm', 'sqp', ...
                           'Display', 'off', ...
                           'MaxIterations', 50, ...
                           'OptimalityTolerance', 1e-3, ...
                           'ConstraintTolerance', 1e-3);
    
    % Solve optimization problem
    cost_function = @(u) simple_cost(u, current_state, d);
    constraint_function = @(u) simple_constraints(u, current_state, obstacles, d);
    
    [u_opt, ~, exitflag] = fmincon(cost_function, x0, [], [], [], [], lb, ub, constraint_function, options);
    
    % Extract first control
    u_next = u_opt(1:d.p.n_u);
    
    % Debug information
    fprintf('Optimization completed with flag %d\n', exitflag);
    fprintf('Control selected: v=%.3f, ω=%.3f\n', u_next(1), u_next(2));
    
    % Handle optimizer failure
    if exitflag <= 0
        warning('Optimizer failed, using fallback control');
        u_next = [0.15; 0.05 * sign(angle_diff)]; % Fallback control
    end
    
catch e
    disp(['Optimization error: ', e.message]);
end
end