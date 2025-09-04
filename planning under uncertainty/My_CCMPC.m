function d = My_CCMPC(obstacles, params)
% d has two components: 1) state: record the state and control inputs and
% time; 2) controller: an object of class optimizer that use ipopt to solve
% optimization in N steps; 3) parameters.

d.p = params;

d.p.x_warm = zeros(d.p.n_x, d.p.N_NMPC+1);
d.p.u_warm = zeros(d.p.n_u, d.p.N_NMPC);

% create NMPC: we need to create an object of class optimizer, this
% requires to save data space for the constraint, objective, optimization
% settings, and controller output
d = create_NMPC(d, obstacles);

% Now, solve the NMPC problem for the simulation time
for t = 1:d.p.sim_steps
    
    [d, u_t, x_t] = solve_NMPC(d,t);

    d = evolve_dynamics(d,t);
    
    % Prepare warm start for next iteration - shift previous solution
    d.p.u_warm(:,1:end-1) = u_t(:,2:end);
    d.p.u_warm(:,end) = u_t(:,end);  % Repeat last control

    d.p.x_warm(:,1) = d.s.x(:,t+1);  % Will be the measured state in next iteration
    d.p.x_warm(:,2:end-1) = x_t(:,3:end);
    d.p.x_warm(:,end) = x_t(:,end);  % Repeat last state

    display(d.s.x(:,t))

end

end

function d = create_NMPC(d, obstacles)
% pre-allocate memory
d.s.x = NaN(d.p.n_x,d.p.sim_steps+1);
d.s.u = NaN(d.p.n_u,d.p.sim_steps);
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
    % control input constraints
    con = [con, d.p.v_min <= u(1,k) <= d.p.v_max];

    % control input constraints
    con = [con, d.p.omega_min <= u(2,k) <= d.p.omega_max];
  
end


% ensure continuity of state trajectory
    % (i.e., direct multiple shooting)
for k = 1:d.p.N_NMPC-1
    if k+1 <= size(x,2)
        con = [con, x(:,k+1) == car_dynamics(x(:,k), u(:,k), d.p.dt)];
    else
        disp(['Error at k=', num2str(k), ': trying to access column ', num2str(k+1), ' of x']);
        break;
    end
end

% add stage cost to objective function
for k = 1:d.p.N_NMPC

%     add stage cost to objective function
    obj = obj + u(:,k)'*d.p.R*u(:,k) + (x(1:2,k)-d.p.goal)'*d.p.Q_k*(x(1:2,k)-d.p.goal);
%     obj = obj + u(:,k)'*d.p.R*u(:,k);
end

% add terminal cost to objective function
obj = obj + (x(1:2, d.p.N_NMPC+1) - d.p.goal)'*d.p.Q*(x(1:2, d.p.N_NMPC+1) - d.p.goal);

%%%%% Now, add the probability constraint %%%%%%%%%
% Error function inverse for probability thresholds
% erf_inv_delta = erfinv(1 - 2*d.p.delta);
% 
% % Uncertainty propagation matrices
% % For simplicity, we use a simplified uncertainty model
% % where the uncertainty grows with time
% Sigma_i = zeros(2, 2, d.p.N_NMPC+1);
% Sigma_i(:,:,1) = [0.01, 0; 0, 0.01]; % Initial position uncertainty
% 
% for k = 2:d.p.N_NMPC+1
%     % Simplified uncertainty propagation (grows with time)
%     Sigma_i(:,:,k) = Sigma_i(:,:,k-1) + diag([0.0005, 0.0005]);
% end

%%%% Compute the probability constraint at stage k among the N horizon 
% for k=2:d.p.N_NMPC+1
%     
%     % for all obstacles;
%     for o = 1:length(obstacles)
%         obstacle = obstacles{o};
% 
%         % Extract robot i's predicted position
%         p_i = x(1:2,k);
%         p_o = obstacle.pos;
% 
%         if obstacle.shape == "circle"
%             % For circular obstacles (simpler case)
%             diff = p_i - p_o;
%             dist = diff' * diff;
% 
%             % Normalized vector
%             a_io = diff / max(dist, 0.01); % Avoid division by zero
% 
%             % Distance
%             b_io = d.p.radius + obstacle.radius;
% 
%             % Combined covariance
%             Sigma_io = Sigma_i(:,:,k) + obstacle.cov;
% 
%             % Probabilistic constraint
%             safety_distance = erf_inv_delta * sqrtm(2*a_io'*Sigma_io*a_io);
%             con = [con, a_io'*(p_i - p_o) - b_io >= safety_distance];
% 
%         elseif obstacle.shape == "ellipse"
%             % For elliptical obstacles (more complex case as per equation 12)
% 
%             % Transformation matrix for the ellipse
%             theta = obstacle.orientation;
%             R_o = [cos(theta), -sin(theta); sin(theta), cos(theta)];
%             S_o = diag(1./obstacle.semi_axes.^2);
%             Omega_io = R_o'*S_o*R_o;
%             Omega_sqrt = sqrtm(Omega_io); % Square root of matrix
% 
%             % Transform to ellipse space
%             p_i_tilde = Omega_sqrt * p_i;
%             p_o_tilde = Omega_sqrt * p_o;
%             diff_tilde = p_i_tilde - p_o_tilde;
%             dist_tilde = diff_tilde'* diff_tilde;
%             
%             % Normalized vector in transformed space
%             a_io_tilde = diff_tilde / max(dist_tilde, 0.01);
% 
%             % In transformed space, the boundary is at distance 1
%             b_io = 1 + d.p.radius*norm(Omega_sqrt*[1;1])/sqrtm(2); % Approximation for transformed radius
% 
%             % Combined covariance in transformed space
%             Sigma_io_tilde = Omega_sqrt * Sigma_i(:,:,k) * Omega_sqrt';
% 
%             % Probabilistic constraint
%             safety_distance = erf_inv_delta * sqrtm(2*a_io_tilde'*Sigma_io_tilde*a_io_tilde);
%             con = [con, a_io_tilde'*(p_i_tilde - p_o_tilde) - b_io >= safety_distance];
%         end
%     end
% end

% Solve optimization problem
ops = sdpsettings('verbose', 0, ...
                  'solver', 'ipopt', ...
                  'usex0', 1);
                  
% Set IPOPT-specific options
ops.ipopt.warm_start_init_point = 'no'; %Primal-only warm start (simpler)
ops.ipopt.mu_strategy = 'adaptive';
ops.ipopt.max_iter = 100;

% % Configure the controller outputs
solutions_out = [u(:); x(:)];

% create MPC controller object
d.c.controller = optimizer(con,obj,ops,{x(:,1), u_warm, x_warm},solutions_out);

d.con = con;
end


function [d, u_t, x_t] = solve_NMPC(d,t)

tic_c = tic;

sol_t = d.c.controller{d.s.x(:,t), d.p.u_warm, d.p.x_warm};

% record CPU time
d.s.CPU_time(t,1) = toc(tic_c);

% extract the control input trajectory
% from the solution of MPC optimization problem
u_t = reshape(sol_t(1:d.p.N_NMPC*2), [d.p.n_u d.p.N_NMPC]);

% extract the state trajectory
% from the solution of MPC optimization problem
x_t = reshape(sol_t(d.p.N_NMPC*2+1:end),[d.p.n_x d.p.N_NMPC+1]);

d.p.x_NMPC_t_1 = x_t;

% assign first element of the solution to the NMPC
% problem as the control input at time t
d.s.u(:,t) = u_t(:,1);
    
end

function x_next = car_dynamics(x, u, dt)
    px = x(1);
    py = x(2);
    psi = x(3);
    
    v = u(1);
    omega = u(2);
    
    % Discretized unicycle dynamics
    x_next(1,:) = px + dt * v * cos(psi*180/pi);
    x_next(2,:) = py + dt* v * sin(psi*180/pi);
    x_next(3,:) = psi + dt * omega;
end

function d = evolve_dynamics(d,t)
    x = d.s.x(:,t);
    u = d.s.u(:,t);
    
    d.s.x(:,t+1) = car_dynamics(x, u, d.p.dt);
end
