function d = test_CCMPC_2(obstacles, params)
d.p = params;

d.p.x_warm = zeros(d.p.n_x, d.p.N_NMPC+1);
d.p.u_warm = zeros(d.p.n_u, d.p.N_NMPC);

% Tolerance for goal reaching (how close car needs to be to goal)
if ~isfield(d.p, 'goal_tolerance')
    d.p.goal_tolerance = 0.1; % Default tolerance, adjust as needed
end

% create NMPC controller
d = create_NMPC(d, obstacles);

% Now, solve the NMPC problem for the simulation time
for t = 1:d.p.sim_steps
    % Check if goal has been reached
    if t > 1
        % Calculate distance to goal
        current_pos = d.s.x(1:2, t);
        distance_to_goal = norm(current_pos - d.p.goal(1:2));
        
        % Calculate orientation difference
        current_orientation = d.s.x(3, t);
        goal_orientation = d.p.goal(3);
        orientation_diff = abs(current_orientation - goal_orientation);

        % Display current distance and orientation difference
        fprintf('Time step %d: Distance to goal = %.4f, Orientation diff = %.4f \n', ...
                t, distance_to_goal, orientation_diff);

        % Check if car is close enough to goal (both position and orientation)
        if distance_to_goal <= d.p.goal_tolerance && orientation_diff <= d.p.orientation_tolerance
            fprintf('Goal reached at time step %d with distance %.4f and orientation diff %.4f!\n', ...
                    t, distance_to_goal, orientation_diff);
            % Set zero control input to stop the car
            d.s.u(:, t:end) = zeros(d.p.n_u, d.p.sim_steps-t+1);
            break;
        end
    end
    
    % Solve NMPC
    [d, u_t, x_t] = solve_NMPC(d, t);

    d = evolve_dynamics(d, t);
    
    % Prepare warm start for next iteration - shift previous solution
    % Check for NaN values and handle them safely
    if any(isnan(u_t(:)))
        disp('Warning: NaN values in control solution at t = ' + string(t));
        % Fill any NaNs with zeros or last valid values
        u_t = fillNaN(u_t, zeros(size(u_t)));
    end
    
    if any(isnan(x_t(:)))
        disp('Warning: NaN values in state solution at t = ' + string(t));
        % Fill any NaNs with appropriate values
        x_t = fillNaN(x_t, repmat(d.s.x(:,t+1), 1, size(x_t,2)));
    end
    
    % Now update the warm start values
    d.p.u_warm(:,1:end-1) = u_t(:,2:end);
    d.p.u_warm(:,end) = u_t(:,end);  % Repeat last control

    d.p.x_warm(:,1) = d.s.x(:,t+1);  % Measured state for next iteration
    
    % Safely copy the predicted states for warm starting
    if size(x_t,2) >= 3
        d.p.x_warm(:,2:end-1) = x_t(:,3:end);
    else
        % Handle the case where x_t doesn't have enough columns
        for i = 2:d.p.N_NMPC
            d.p.x_warm(:,i) = d.s.x(:,t+1); % Use current state as fallback
        end
    end
    
    d.p.x_warm(:,end) = x_t(:,end);  % Repeat last state
    
    % Verify warm start values are valid
    if any(isnan(d.p.u_warm(:))) || any(isnan(d.p.x_warm(:)))
        disp('Error: Warm start values contain NaN after preparation');
        % Set all warm start values to safe defaults
        d.p.u_warm = zeros(d.p.n_u, d.p.N_NMPC);
        d.p.x_warm = repmat(d.s.x(:,t+1), 1, d.p.N_NMPC+1);
    end

    display(d.s.x(:,t))
end

% Display final state
if exist('t', 'var')
    disp('Final state:');
    display(d.s.x(:,t));
end
end

% Helper function to fill NaN values
function data_filled = fillNaN(data, replacement)
    data_filled = data;
    nan_indices = isnan(data);
    data_filled(nan_indices) = replacement(nan_indices);
end

function d = create_NMPC(d, obstacles)
% pre-allocate memory
d.s.x = zeros(d.p.n_x, d.p.sim_steps+1);
d.s.u = zeros(d.p.n_u, d.p.sim_steps);
d.s.CPU_time = NaN(d.p.n_u, d.p.sim_steps);

% set initial state
d.s.x(:,1) = d.p.x0;

% Store problem data in the controller structure
d.c.obstacles = obstacles;
d.c.erf_inv_delta = erfinv(1 - 2*d.p.delta);  % Error function inverse for probability
d.c.n_x = d.p.n_x;
d.c.n_u = d.p.n_u;
d.c.N_NMPC = d.p.N_NMPC;
d.c.dt = d.p.dt;
d.c.R = d.p.R;
d.c.Q = d.p.Q;
d.c.Q_k = d.p.Q_k;
d.c.goal = d.p.goal;
d.c.radius = d.p.radius;
d.c.Sigma_0 = d.p.Sigma_0;
d.c.Sigma_propogagtion = d.p.Sigma_propogagtion;
d.c.v_max = d.p.v_max;
d.c.omega_max = d.p.omega_max;

end

function [d, u_t, x_t] = solve_NMPC(d, t)
tic_c = tic;

try
    % Initialize the optimization variables
    n_vars = d.p.n_u * d.p.N_NMPC + d.p.n_x * (d.p.N_NMPC + 1);
    
    % Use the warm start as the initial guess for optimization
    x0 = [reshape(d.p.u_warm, [], 1); reshape(d.p.x_warm, [], 1)];
    
    % Define lower and upper bounds for the variables
    % Control bounds
    lb_u = repmat([-d.p.v_max; -d.p.omega_max], d.p.N_NMPC, 1);
    ub_u = repmat([d.p.v_max; d.p.omega_max], d.p.N_NMPC, 1);
    
    % State bounds (can be set to large values if needed)
    lb_x = repmat([-inf; -inf; -inf], d.p.N_NMPC + 1, 1);
    ub_x = repmat([inf; inf; inf], d.p.N_NMPC + 1, 1);
    
    % Combine bounds
    lb = [lb_u; lb_x];
    ub = [ub_u; ub_x];
    
    % Initial state constraint
    Aeq = zeros(d.p.n_x, n_vars);
    Aeq(:, d.p.n_u * d.p.N_NMPC + (1:d.p.n_x)) = eye(d.p.n_x);
    beq = d.s.x(:, t);
    
    % Set optimization options
    options = optimoptions('fmincon', ...
                           'Algorithm', 'interior-point', ...
                           'Display', 'off', ...
                           'MaxIterations', 100, ...
                           'SpecifyObjectiveGradient', false, ...
                           'SpecifyConstraintGradient', false, ...
                           'OptimalityTolerance', 1e-4, ...
                           'ConstraintTolerance', 1e-4);
    
    % Solve the optimization problem using fmincon
    [sol_t, ~, exitflag] = fmincon(@(z) objective_function(z, d, d.c.obstacles, t), x0, [], [], Aeq, beq, lb, ub, ...
                                   @(z) nonlinear_constraints(z, d, d.c.obstacles, t), options);
    
    % Check if optimization was successful
    if exitflag <= 0
        warning('Optimization problem status: %d', exitflag);
    end
    
    % Record CPU time
    d.s.CPU_time(t,1) = toc(tic_c);
    
    % Extract the control input trajectory from the solution
    u_t = reshape(sol_t(1:d.p.n_u*d.p.N_NMPC), [d.p.n_u, d.p.N_NMPC]);
    
    % Extract the state trajectory from the solution
    x_t = reshape(sol_t(d.p.n_u*d.p.N_NMPC+1:end), [d.p.n_x, d.p.N_NMPC+1]);
    
    % Store predicted trajectory
    d.p.x_NMPC_t_1 = x_t;
    
    % Assign first element of the solution as the control input at time t
    d.s.u(:,t) = u_t(:,1);
    
catch e
    disp(['Error in solve_NMPC at t = ', num2str(t)]);
    disp(e.message);
    
    % Provide fallback values in case of error
    if t > 1
        % Use previous control as fallback
        d.s.u(:,t) = d.s.u(:,t-1);
        u_t = repmat(d.s.u(:,t), 1, d.p.N_NMPC);
        x_t = repmat(d.s.x(:,t), 1, d.p.N_NMPC+1);
    else
        % For t=1, use zeros as fallback
        d.s.u(:,t) = zeros(d.p.n_u, 1);
        u_t = zeros(d.p.n_u, d.p.N_NMPC);
        x_t = repmat(d.s.x(:,t), 1, d.p.N_NMPC+1);
    end
    
    d.s.CPU_time(t,1) = toc(tic_c);
end
    
end

function [cost] = objective_function(z, d, obstacles, t)
    % Extract control inputs and states from z
    n_u = d.p.n_u;
    n_x = d.p.n_x;
    N = d.p.N_NMPC;
    
    u = reshape(z(1:n_u*N), [n_u, N]);
    x = reshape(z(n_u*N+1:end), [n_x, N+1]);
    
    % Initialize cost
    cost = 0;
    
    % Stage cost
    for k = 1:N
        % Control cost
        cost = cost + u(:,k)'*d.p.R*u(:,k);
        
        % State cost (tracking to goal)
        cost = cost + (x(:,k) - d.p.goal)'*d.p.Q_k*(x(:,k) - d.p.goal);
        
        % Get predicted distance to goal
        predicted_pos = x(1:2,k);
        distance_to_goal = norm(predicted_pos - d.p.goal(1:2));
        
        % Penalty for velocity near goal (soft constraint)
        if k > 1
            scaling_factor = min(distance_to_goal / 1.0, 1.0);
            max_vel_near_goal = d.p.v_max * scaling_factor;
            max_omega_near_goal = d.p.omega_max * scaling_factor;
            
            % Soft constraints as penalties
            if u(1,k) > max_vel_near_goal
                cost = cost + 1000 * (u(1,k) - max_vel_near_goal)^2;
            end
            
            if u(2,k) > max_omega_near_goal
                cost = cost + 1000 * (u(2,k) - max_omega_near_goal)^2;
            end
        end
        
        % Collision avoidance constraints are now handled in nonlinear_constraints
        % We removed the soft penalties from the objective function
    end
    
    % Terminal cost
    cost = cost + (x(:,N+1) - d.p.goal)'*d.p.Q*(x(:,N+1) - d.p.goal);
end

function [c, ceq] = nonlinear_constraints(z, d, obstacles, t)
    % Extract control inputs and states from z
    n_u = d.p.n_u;
    n_x = d.p.n_x;
    N = d.p.N_NMPC;
    
    u = reshape(z(1:n_u*N), [n_u, N]);
    x = reshape(z(n_u*N+1:end), [n_x, N+1]);
    
    % System dynamics constraints (equality constraints)
    ceq = zeros(n_x * N, 1);
    for k = 1:N
        next_state = car_dynamics(x(:,k), u(:,k), d.p.dt);
        ceq((k-1)*n_x+1:k*n_x) = x(:,k+1) - next_state;
    end
    
    % Collision avoidance constraints (inequality constraints)
    % Count the number of obstacles and prediction steps
    n_obstacles = length(obstacles);
    n_constraints = 0;
    for k = 2:N+1  % Start from 2 to skip initial state
        n_constraints = n_constraints + n_obstacles;
    end
    
    % Allocate space for inequality constraints
    c = zeros(n_constraints, 1);
    
    % Fill in constraints
    constraint_idx = 1;
    for k = 2:N+1  % Start from 2 to skip initial state
        % Compute uncertainty for this prediction step
        Sigma_i = zeros(2, 2);
        if k == 1
            Sigma_i = d.p.Sigma_0;
        else
            Sigma_i = d.p.Sigma_0 + (k-1)*d.p.Sigma_propogagtion;
        end
        
        % Check constraints for each obstacle
        for o = 1:n_obstacles
            obstacle = obstacles{o};
            
            % Extract robot position
            p_i = x(1:2,k);
            p_o = obstacle.pos;
            
            % Relative position
            p_io = p_i - p_o;
            
            % Normal 
            a_io = p_io / norm(p_io);
            
             if obstacle.shape == "circle"
                % Distance
                b_io = d.p.radius + obstacle.radius;

                % Combined covariance
                Sigma_io = Sigma_i + obstacle.cov;

                % Probabilistic constraint
                quadratic_term = 2*a_io'*Sigma_io*a_io;
                safety_distance = d.c.erf_inv_delta *sqrt(max(quadratic_term, 1e-6));  % Avoid numerical issues
                
                % Constraint: a_io' * p_io - b_io >= safety_distance
                % Written as: safety_distance - (a_io' * p_io - b_io) <= 0

                c(constraint_idx) = safety_distance - ( a_io' * p_io - b_io);
                constraint_idx = constraint_idx + 1;

            elseif obstacle.shape == "ellipse"
                % For elliptical obstacles (more complex case as per equation 12)

                % Transformation matrix for the ellipse
                theta = obstacle.orientation;
                R_o = [cos(theta), -sin(theta); sin(theta), cos(theta)];

                S_io = diag(1./(obstacle.semi_axes + d.p.radius).^2);
                Omega_io = R_o'*S_io*R_o;
                Omega_sqrt = sqrtm(Omega_io); % Square root of matrix

                % In transformed space, the boundary is at distance 1
                b_io = 1;

                % Combined covariance in transformed space
                Sigma_io_tilde = Omega_sqrt * (Sigma_i +R_o' * obstacle.cov * R_o) * Omega_sqrt';
                
                % Probabilistic constraint
                quadratic_term = 2*a_io'*Sigma_io_tilde*a_io;
                safety_distance = d.c.erf_inv_delta *sqrt(max(quadratic_term, 1e-6));  % Avoid numerical issues

                c(constraint_idx) = safety_distance - ( a_io' * Omega_sqrt * p_io - b_io);
                constraint_idx = constraint_idx + 1;
             end
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
    x_next(1,:) = px + dt * v * cos(psi);
    x_next(2,:) = py + dt * v * sin(psi);
    x_next(3,:) = psi + dt * omega;
end

function d = evolve_dynamics(d, t)
    x = d.s.x(:,t);
    u = d.s.u(:,t);
    
    d.s.x(:,t+1) = car_dynamics(x, u, d.p.dt);
end