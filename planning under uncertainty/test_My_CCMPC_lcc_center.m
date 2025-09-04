function d = test_My_CCMPC_lcc_center(obstacles, params)
d.p = params;

d.p.x_warm = zeros(d.p.n_x, d.p.N_NMPC+1);
d.p.u_warm = zeros(d.p.n_u, d.p.N_NMPC);

% Tolerance for goal reaching (how close car needs to be to goal)
if ~isfield(d.p, 'goal_tolerance')
    d.p.goal_tolerance = 0.1; % Default tolerance, adjust as needed
end

% create NMPC: we need to create an object of class optimizer, this
% requires to save data space for the constraint, objective, optimization
% settings, and controller output
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
    
    % Solve NMPC with proper handling of output parameters
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
d.s.x = zeros(d.p.n_x,d.p.sim_steps+1);
d.s.u = zeros(d.p.n_u,d.p.sim_steps);
d.s.CPU_time = NaN(d.p.n_u,d.p.sim_steps);

% In your simulation loop
x_warm = sdpvar(d.p.n_x, d.p.N_NMPC+1); % Initial warm start for states
u_warm = sdpvar(d.p.n_u, d.p.N_NMPC);   % Initial warm start for control inputs

% set initial state
% (d.p.x0 is an argument of the main function)
d.s.x(:,1) = d.p.x0;

% define control inputs as decision variables
u = sdpvar(d.p.n_u,d.p.N_NMPC,'full');

% define states as decision variables
x = sdpvar(d.p.n_x,d.p.N_NMPC + 1,'full');

% initialize objective function
obj = 0;
    
% initialize constraints
con = [];

% define constraints for the control variable
for k = 1:d.p.N_NMPC
    
    % Calculate predicted distance to goal for this state
    predicted_pos = x(1:2,k);
    distance_to_goal = (predicted_pos - d.p.goal(1:2))' * (predicted_pos - d.p.goal(1:2));

    % Add velocity reduction constraint near goal
    if k > 1  % Skip first step for feasibility
        scaling_factor = min(distance_to_goal / 1.0, 1.0);
        max_vel_near_goal = d.p.v_max * scaling_factor;
        % Add a soft constraint for velocity reduction (to ensure feasibility)
        slack = sdpvar(1);
        con = [con, u(1,k) <= max_vel_near_goal + slack, slack >= 0];
        max_omega_near_goal = d.p.omega_max * scaling_factor;
        con = [con, u(2,k) <= max_omega_near_goal + slack, slack >= 0];
        obj = obj + 1000*slack; % Penalize violation of this soft constraint
    end
      
end

% ensure continuity of state trajectory
% (i.e., direct multiple shooting)
for k = 1:d.p.N_NMPC
    con = [con, x(:,k+1) == car_dynamics(x(:,k), u(:,k), d.p.dt)];
end

% add probabilistic collision probability constraint
% Error function inverse for probability thresholds
erf_inv_delta = erfinv(1 - 2*d.p.delta);

% Uncertainty propagation matrices
% For simplicity, we use a simplified uncertainty model
% where the uncertainty grows with time
Sigma_i = zeros(2, 2, d.p.N_NMPC+1);
Sigma_i(:,:,1) = d.p.Sigma_0;

for k = 2:d.p.N_NMPC+1
    % Simplified uncertainty propagation (grows with time)
    Sigma_i(:,:,k) = Sigma_i(:,:,k-1) + d.p.Sigma_propogagtion;
end
% 
%%% Compute the probability constraint at stage k among the N horizon 
for k=2:d.p.N_NMPC+1
    
    % for all obstacles;
    for o = 1:length(obstacles)
        obstacle = obstacles{o};

        % Extract robot i's predicted position
        p_i = x(1:2,k);
        p_o = obstacle.pos;

        % For circular obstacles (simpler case)
        p_io = p_i - p_o;
        
        % Distance
        b_io = d.p.radius + obstacle.radius;

        % Combined covariance
        Sigma_io = Sigma_i(:,:,k) + obstacle.cov;

        % Probabilistic constraint
        quadratic_term = 2*p_io'*Sigma_io*p_io;
        safety_distance = erf_inv_delta * (quadratic_term^0.5);
    
        slack_collision = sdpvar(1);
        con = [con,  p_io' * p_io - b_io * sqrtm(p_io' * p_io) >= safety_distance - slack_collision, slack_collision >= 0];
        obj = obj + 5000*slack_collision; % High penalty for violation
       
    end
end

% add stage cost to objective function
for k = 1:d.p.N_NMPC
       
    % add stage cost to objective function
    obj = obj + u(:,k)'*d.p.R*u(:,k) + (x(:, d.p.N_NMPC+1) - d.p.goal)'*d.p.Q_k*(x(:, d.p.N_NMPC+1) - d.p.goal);
end

% add terminal cost to objective function
obj = obj + (x(:, d.p.N_NMPC+1) - d.p.goal)'*d.p.Q*(x(:, d.p.N_NMPC+1) - d.p.goal);

% Solve optimization problem
ops = sdpsettings('verbose', 0, ...
                  'solver', 'ipopt', ...
                  'usex0', 1);
                  
% Set IPOPT-specific options
ops.ipopt.warm_start_init_point = 'no'; % Primal-only warm start (simpler)
ops.ipopt.mu_strategy = 'adaptive';
ops.ipopt.max_iter = 100;


% % Configure the controller outputs
solutions_out = [u(:); x(:)];

% Create MPC controller object - use additional outputs to fix warning
d.c.controller = optimizer(con, obj, ops, {x(:,1), u_warm, x_warm}, solutions_out);

d.con = con;
end

function [d, u_t, x_t] = solve_NMPC(d, t)
tic_c = tic;

try
    % Use the full output syntax to avoid warnings about initial guesses
    [sol_t, problem, ~, ~, P] = d.c.controller{d.s.x(:,t), d.p.u_warm, d.p.x_warm};
    
    % Store the optimizer object P for future warm starts
    if ~isempty(P)
        d.c.P = P;
    end
    
    % Check if optimization was successful
    if problem ~= 0
        error('Optimization problem status: %d', problem);
    end
    
    % record CPU time
    d.s.CPU_time(t,1) = toc(tic_c);
    
    % extract the control input trajectory
    % from the solution of MPC optimization problem
    u_t = reshape(sol_t(1:d.p.N_NMPC*2), [d.p.n_u d.p.N_NMPC]);
    
    % extract the state trajectory
    % from the solution of MPC optimization problem
    x_t = reshape(sol_t(d.p.N_NMPC*2+1:end), [d.p.n_x d.p.N_NMPC+1]);
    
    d.p.x_NMPC_t_1 = x_t;
    
    % assign first element of the solution to the NMPC
    % problem as the control input at time t
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

function x_next = car_dynamics(x, u, dt)
    px = x(1);
    py = x(2);
    psi = x(3);
    
    v = u(1);
    omega = u(2);
    
    % Discretized unicycle dynamics
    x_next(1,:) = px + dt * v * cos(psi);
    x_next(2,:) = py + dt* v * sin(psi);
    x_next(3,:) = psi + dt * omega;
end

function d = evolve_dynamics(d,t)
    x = d.s.x(:,t);
    u = d.s.u(:,t);
    
    d.s.x(:,t+1) = car_dynamics(x, u, d.p.dt);
end