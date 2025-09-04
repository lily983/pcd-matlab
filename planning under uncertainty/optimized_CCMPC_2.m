function d = optimized_CCMPC_2(obstacles, params)
d.p = params;

% Initialize warm start variables
d.p.x_warm = zeros(d.p.n_x, d.p.N_NMPC+1);
d.p.u_warm = zeros(d.p.n_u, d.p.N_NMPC);
d.c.obstacles = obstacles;  % Store obstacles data

% Set erf inverse delta for better performance
d.c.erf_inv_delta = sqrt(2) * erfinv(1 - 2*d.p.delta);

% Pre-compute obstacle data to avoid repetitive calculations
for o = 1:length(obstacles)
    obstacle = obstacles{o};
    theta = obstacle.orientation;
    R_o = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    S_io = diag(1./(obstacle.semi_axes + d.p.radius).^2);
    obstacles{o}.R_o = R_o;
    obstacles{o}.S_io = S_io;
    obstacles{o}.Omega_io = R_o'*S_io*R_o;
    obstacles{o}.Omega_sqrt = sqrtm(obstacles{o}.Omega_io);
end
d.c.obstacles = obstacles;

% Set up goal reaching tolerance
if ~isfield(d.p, 'goal_tolerance')
    d.p.goal_tolerance = 0.1;
end

% create NMPC controller
d = create_NMPC(d);

% Set up obstacle repulsion field (new)
d.p.obstacle_repulsion = 3.0;  % Strength of the repulsion field
d.p.repulsion_radius = 1.0;    % Radius of influence for repulsion
d.p.safe_distance = 0.5;       % Distance considered "safe" from obstacles

% Now, solve the NMPC problem for the simulation time
for t = 1:d.p.sim_steps
    
    % Add to solve_NMPC function, before the try block
    if t == 1
        % Force non-zero initial velocity to get things started
        d.p.u_warm(1,:) = ones(1,d.p.N_NMPC_current) * 0.1; 
    end
    % Check if goal has been reached
    if t > 1
        % Calculate distance to goal
        current_pos = d.s.x(1:2, t);
        distance_to_goal = norm(current_pos - d.p.goal(1:2));
        
        % Calculate orientation difference
        current_orientation = d.s.x(3, t);
        goal_orientation = d.p.goal(3);
        orientation_diff = abs(current_orientation - goal_orientation);
        orientation_diff = min(orientation_diff, 2*pi - orientation_diff);

        % Display current distance and orientation difference
        fprintf('Time step %d: Distance to goal = %.4f, Orientation diff = %.4f, CPU time = %.4f\n', ...
                t, distance_to_goal, orientation_diff, d.s.CPU_time(t-1));

        % Check if car is close enough to goal (both position and orientation)
        if distance_to_goal <= d.p.goal_tolerance && orientation_diff <= d.p.orientation_tolerance
            fprintf('Goal reached at time step %d with distance %.4f and orientation diff %.4f!\n', ...
                    t, distance_to_goal, orientation_diff);
            % Set zero control input to stop the car
            d.s.u(:, t:end) = zeros(d.p.n_u, d.p.sim_steps-t+1);
            break;
        end
    end
    
    % If we're near obstacles, increase the prediction horizon adaptively
    if t > 1
        obstacle_distance = min_obstacle_distance(d.s.x(:,t), obstacles);
        if obstacle_distance < d.p.repulsion_radius && d.p.adaptive_horizon
            % Increase horizon near obstacles
            d.p.N_NMPC_current = min(d.p.N_NMPC_max, d.p.N_NMPC_base + round((d.p.repulsion_radius - obstacle_distance) * 10));
            fprintf('Increased horizon to %d at distance %.2f\n', d.p.N_NMPC_current, obstacle_distance);
        else
            d.p.N_NMPC_current = d.p.N_NMPC_base;
        end
    else
        d.p.N_NMPC_current = d.p.N_NMPC_base;
    end
    
    % Solve NMPC
    [d, u_t, x_t] = solve_NMPC(d, t);

    % Evolve dynamics
    d = evolve_dynamics(d, t);
    
    % Prepare warm start for next iteration
    d = prepare_warm_start(d, t, u_t, x_t);
end

% Display final state
if exist('t', 'var')
    disp('Final state:');
    disp(d.s.x(:,t));
    
    % Display performance metrics
    avg_cpu_time = mean(d.s.CPU_time(1:t));
    max_cpu_time = max(d.s.CPU_time(1:t));
    fprintf('Average CPU time: %.4f seconds\n', avg_cpu_time);
    fprintf('Maximum CPU time: %.4f seconds\n', max_cpu_time);
end
end

function d = create_NMPC(d)
% pre-allocate memory
d.s.x = zeros(d.p.n_x, d.p.sim_steps+1);
d.s.u = zeros(d.p.n_u, d.p.sim_steps);
d.s.CPU_time = zeros(d.p.sim_steps, 1);
d.s.tangent_help = zeros(3,1);
d.s.obstacle_distances = zeros(d.p.sim_steps, 1);

% set initial state
d.s.x(:,1) = d.p.x0;
end

function [d, u_t, x_t] = solve_NMPC(d, t)
tic_c = tic;

try
    % Get the current prediction horizon (may be adaptive)
    N = d.p.N_NMPC_current;
    
    % Initialize the optimization variables
    n_vars = d.p.n_u * N + d.p.n_x * (N + 1);
    
    % Resize warm start if needed due to adaptive horizon
    if size(d.p.u_warm, 2) ~= N
        % Adjust warm start to match current horizon
        old_horizon = size(d.p.u_warm, 2);
        if old_horizon < N
            d.p.u_warm = [d.p.u_warm, repmat(d.p.u_warm(:,end), 1, N-old_horizon)];
            d.p.x_warm = [d.p.x_warm, repmat(d.p.x_warm(:,end), 1, N-old_horizon+1)];
        else
            d.p.u_warm = d.p.u_warm(:,1:N);
            d.p.x_warm = d.p.x_warm(:,1:N+1);
        end
    end
    
    % Use the warm start as the initial guess for optimization
    x0 = [reshape(d.p.u_warm, [], 1); reshape(d.p.x_warm, [], 1)];
    
    % Define lower and upper bounds for the variables
    lb_u = repmat([d.p.v_min; d.p.omega_min], N, 1);
    ub_u = repmat([d.p.v_max; d.p.omega_max], N, 1);
    
    % State bounds 
    lb_x = repmat([-inf; -inf; -inf], N + 1, 1);
    ub_x = repmat([inf; inf; inf], N + 1, 1);
    
    % Combine bounds
    lb = [lb_u; lb_x];
    ub = [ub_u; ub_x];
    
    % Initial state constraint
    Aeq = zeros(d.p.n_x, n_vars);
    Aeq(:, d.p.n_u * N + (1:d.p.n_x)) = eye(d.p.n_x);
    beq = d.s.x(:, t);
    
    % Current distance to obstacles for adaptive parameters
    current_obstacle_distance = min_obstacle_distance(d.s.x(:,t), d.c.obstacles);
    d.s.obstacle_distances(t) = current_obstacle_distance;
    
    % Adaptive optimization settings based on obstacle proximity
    if current_obstacle_distance < d.p.safe_distance
        % Near obstacles: tighter tolerances, more iterations
        options = optimoptions('fmincon', ...
                           'Algorithm', 'sqp', ... % SQP often better with nonlinear constraints
                           'Display', 'off', ...
                           'MaxIterations', 200, ... % More iterations
                           'MaxFunctionEvaluations', 5000, ... % More function evaluations
                           'SpecifyObjectiveGradient', false, ...
                           'SpecifyConstraintGradient', false, ...
                           'OptimalityTolerance', 1e-5, ... % Tighter tolerances
                           'ConstraintTolerance', 1e-5);
    else
        % Away from obstacles: looser tolerances for speed
        options = optimoptions('fmincon', ...
                           'Algorithm', 'interior-point', ...
                           'Display', 'off', ...
                           'MaxIterations', 100, ...
                           'MaxFunctionEvaluations', 3000, ...
                           'SpecifyObjectiveGradient', false, ...
                           'SpecifyConstraintGradient', false, ...
                           'OptimalityTolerance', 1e-4, ...
                           'ConstraintTolerance', 1e-4);
    end
    
    % Solve the optimization problem using fmincon
    [sol_t, ~, exitflag, output] = fmincon(@(z) objective_function(z, d, t, current_obstacle_distance), ...
                                   x0, [], [], Aeq, beq, lb, ub, ...
                                   @(z) nonlinear_constraints(z, d, t), options);
    
    % Check if optimization was successful
    if exitflag <= 0
        warning('Optimization problem status: %d at time %d, iterations: %d', exitflag, t, output.iterations);
        % If failure and we're not already at minimum horizon, try again with shorter horizon
        if N > 5 && d.p.adaptive_horizon
            d.p.N_NMPC_current = max(5, N - 2);
            fprintf('Reducing horizon to %d and retrying\n', d.p.N_NMPC_current);
            [d, u_t, x_t] = solve_NMPC(d, t);
            return;
        end
    end
    
    % Record CPU time
    d.s.CPU_time(t) = toc(tic_c);
    
    % Extract the control input trajectory from the solution
    u_t = reshape(sol_t(1:d.p.n_u*N), [d.p.n_u, N]);
    
    % Extract the state trajectory from the solution
    x_t = reshape(sol_t(d.p.n_u*N+1:end), [d.p.n_x, N+1]);
    
    % Store predicted trajectory
    d.p.x_NMPC_t_1 = x_t;
    
    % Assign first element of the solution as the control input at time t
    d.s.u(:,t) = u_t(:,1);
    
catch e
    disp(['Error in solve_NMPC at t = ', num2str(t)]);
    disp(e.message);
    
    % Provide fallback values in case of error
    if t > 1
        % Use previous control as fallback but with reduced velocity
        d.s.u(:,t) = d.s.u(:,t-1) * 0.7; % Slow down
        u_t = repmat(d.s.u(:,t), 1, d.p.N_NMPC_current);
        x_t = repmat(d.s.x(:,t), 1, d.p.N_NMPC_current+1);
    else
        % For t=1, use small forward velocity as fallback
        d.s.u(:,t) = [0.1; 0]; % Small forward velocity
        u_t = zeros(d.p.n_u, d.p.N_NMPC_current);
        x_t = repmat(d.s.x(:,t), 1, d.p.N_NMPC_current+1);
    end
    
    d.s.CPU_time(t) = toc(tic_c);
end
    
end

function [cost] = objective_function(z, d, t, current_obstacle_distance)
    % Extract control inputs and states from z
    n_u = d.p.n_u;
    n_x = d.p.n_x;
    N = d.p.N_NMPC_current;
    
    u = reshape(z(1:n_u*N), [n_u, N]);
    x = reshape(z(n_u*N+1:end), [n_x, N+1]);
    
    % Initialize cost
    cost = 0;
    
    % Adaptive weights based on distance to obstacles
    obstacle_factor = max(0, 1 - current_obstacle_distance/d.p.repulsion_radius);
    Q_obstacle = d.p.Q_k;
    
    % If very close to obstacle, prioritize obstacle avoidance over goal reaching
    if current_obstacle_distance < d.p.safe_distance
        % Reduce weight on goal tracking when close to obstacles
        goal_weight_factor = max(0.2, current_obstacle_distance/d.p.safe_distance);
        Q_goal = d.p.Q_k * goal_weight_factor;
    else
        Q_goal = d.p.Q_k;
    end
    
    % First-derivative penalty on control to encourage smooth changes
    u_prev = [zeros(n_u, 1), u(:,1:end-1)]; 
    if t > 1
        u_prev(:,1) = d.s.u(:,t-1);
    end
    
    % Stage cost
    for k = 1:N
        % Control cost
        cost = cost + u(:,k)'*d.p.R*u(:,k);
        
        % Control rate cost (smoothness)
        control_rate = u(:,k) - u_prev(:,k);
        cost = cost + control_rate'*diag([0.2, 0.2])*control_rate;
        
        % State cost (tracking to goal)
        cost = cost + (x(:,k) - d.p.goal)'*Q_goal*(x(:,k) - d.p.goal);
        
        % Obstacle avoidance soft cost (new)
        for o = 1:length(d.c.obstacles)
            obstacle = d.c.obstacles{o};
            p_i = x(1:2,k);
            p_o = obstacle.pos;
            
            % Distance to obstacle center
            dist_to_obs = norm(p_i - p_o) - norm(obstacle.semi_axes);
            
            % Apply repulsive potential when within influence radius
            if dist_to_obs < d.p.repulsion_radius
                repulsion = d.p.obstacle_repulsion * (1 - dist_to_obs/d.p.repulsion_radius)^2;
                cost = cost + repulsion;
            end
        end
        
        % Get predicted distance to goal
        predicted_pos = x(1:2,k);
        distance_to_goal = norm(predicted_pos - d.p.goal(1:2));
        
        % Adaptive velocity control near goal
        if distance_to_goal < 0.5
            % Progressive velocity reduction as we approach the goal
            scaling_factor = max(0.1, distance_to_goal / 0.5);
            target_velocity = d.p.v_max * scaling_factor;
            
            % Penalize excess velocity near goal
            if u(1,k) > target_velocity
                cost = cost + 10 * (u(1,k) - target_velocity)^2;
            end
        end
    end
    
    % Terminal cost
    cost = cost + (x(:,N+1) - d.p.goal)'*d.p.Q*(x(:,N+1) - d.p.goal);
end

function [c, ceq] = nonlinear_constraints(z, d, t)
    % Extract control inputs and states from z
    n_u = d.p.n_u;
    n_x = d.p.n_x;
    N = d.p.N_NMPC_current;
    
    u = reshape(z(1:n_u*N), [n_u, N]);
    x = reshape(z(n_u*N+1:end), [n_x, N+1]);
    
    % System dynamics constraints (equality constraints)
    ceq = zeros(n_x * N, 1);
    for k = 1:N
        next_state = car_dynamics(x(:,k), u(:,k), d.p.dt);
        ceq((k-1)*n_x+1:k*n_x) = x(:,k+1) - next_state;
    end
    
    % Collision avoidance constraints (inequality constraints)
    obstacles = d.c.obstacles;
    n_obstacles = length(obstacles);
    
    % Skip first few steps for constraint checking to improve performance
    k_start = 1; % Can be increased to skip initial steps if needed
    
    % Calculate number of constraints
    n_constraints = n_obstacles * (N+1-k_start);
    
    % Allocate space for inequality constraints
    c = zeros(n_constraints, 1);
    
    % Fill in constraints
    constraint_idx = 1;
    for k = k_start:N+1
        % Compute uncertainty for this prediction step
        if k == 1
            Sigma_i = d.p.Sigma_0;
        else
            Sigma_i = d.p.Sigma_0 + (k-1)*d.p.Sigma_propogagtion;
        end
        
        for o = 1:n_obstacles
            obstacle = obstacles{o};
            
            % Extract robot position and obstacle info
            p_i = x(1:2,k);
            p_o = obstacle.pos;
            
            % Relative position
            p_io = p_i - p_o;
            
            % Skip constraint if very far from obstacle (optimization)
            if norm(p_io) > 2*(norm(obstacle.semi_axes) + d.p.radius + 0.5)
                c(constraint_idx) = -1; % Automatically satisfied
                constraint_idx = constraint_idx + 1;
                continue;
            end
            
            % Use pre-computed obstacle matrices
            R_o = obstacle.R_o;
            Omega_sqrt = obstacle.Omega_sqrt;
            
            % Normal vector (for stability, ensure it's normalized)
            p_io_norm = norm(p_io);
            if p_io_norm < 1e-6
                % If very close to obstacle center, move slightly away
                p_io = [1e-6; 1e-6];
                p_io_norm = norm(p_io);
            end
            a_io = p_io / p_io_norm;
            
            % Combined covariance
            Sigma_io = Sigma_i + obstacle.cov;
            
            % Transformed covariance
            Sigma_io_tilde = Omega_sqrt * Sigma_io * Omega_sqrt';
            
            % Probabilistic constraint (simplified calculation)
            quadratic_term = 2*a_io'*Sigma_io_tilde*a_io;
            transformed_dist = a_io' * Omega_sqrt * p_io - 1;
            
            % Fast approximation for collision probability
            prob = 0.5 - 0.5*erf(transformed_dist / sqrt(max(quadratic_term, 1e-6)));
            
            % Apply constraint
            c(constraint_idx) = prob - d.p.delta;
            constraint_idx = constraint_idx + 1;
        end
    end
end

function x_next = car_dynamics(x, u, dt)
    px = x(1);
    py = x(2);
    psi = x(3);
    
    v = u(1);
    omega = u(2);
    
    % Discretized unicycle dynamics
    x_next = zeros(size(x));
    x_next(1) = px + dt * v * cos(psi);
    x_next(2) = py + dt * v * sin(psi);
    x_next(3) = psi + dt * omega;
end

function d = evolve_dynamics(d, t)
    x = d.s.x(:,t);
    u = d.s.u(:,t);
    
    d.s.x(:,t+1) = car_dynamics(x, u, d.p.dt);
end

function d = prepare_warm_start(d, t, u_t, x_t)
    % Check for NaN values and handle them safely
    if any(isnan(u_t(:)))
        disp('Warning: NaN values in control solution at t = ' + string(t));
        u_t = fillNaN(u_t, zeros(size(u_t)));
    end
    
    if any(isnan(x_t(:)))
        disp('Warning: NaN values in state solution at t = ' + string(t));
        x_t = fillNaN(x_t, repmat(d.s.x(:,t+1), 1, size(x_t,2)));
    end
    
    % Now update the warm start values with proper shifting
    N = d.p.N_NMPC_current;
    
    % Shift controls
    if size(u_t, 2) > 1
        d.p.u_warm(:,1:N-1) = u_t(:,2:end);
        d.p.u_warm(:,N) = u_t(:,end);  % Repeat last control
    else
        d.p.u_warm = repmat(u_t, 1, N);
    end

    % Update states
    d.p.x_warm(:,1) = d.s.x(:,t+1);  % Measured state for next iteration
    
    if size(x_t, 2) > 2
        d.p.x_warm(:,2:N) = x_t(:,3:end);
        d.p.x_warm(:,N+1) = x_t(:,end);  % Repeat last state
    else
        % Not enough predicted states, initialize with current state
        d.p.x_warm = repmat(d.s.x(:,t+1), 1, N+1);
    end
    
    % Verify warm start values are valid
    if any(isnan(d.p.u_warm(:))) || any(isnan(d.p.x_warm(:)))
        disp('Error: Warm start values contain NaN after preparation');
        % Set all warm start values to safe defaults
        d.p.u_warm = zeros(d.p.n_u, N);
        d.p.x_warm = repmat(d.s.x(:,t+1), 1, N+1);
    end
end

% Helper function to fill NaN values
function data_filled = fillNaN(data, replacement)
    data_filled = data;
    nan_indices = isnan(data);
    data_filled(nan_indices) = replacement(nan_indices);
end

% Function to calculate minimum distance to any obstacle
function min_dist = min_obstacle_distance(state, obstacles)
    px = state(1);
    py = state(2);
    pos = [px; py];
    
    min_dist = inf;
    for i = 1:length(obstacles)
        obstacle = obstacles{i};
        p_o = obstacle.pos;
        
        % Basic Euclidean distance to obstacle center minus obstacle size
        % This is a rough approximation but works for quick calculations
        dist = norm(pos - p_o) - max(obstacle.semi_axes);
        min_dist = min(min_dist, dist);
    end
end