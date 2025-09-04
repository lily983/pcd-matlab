%% Probabilistic Collision Avoidance for 2D Car Model
% Based on "Chance-Constrained Collision Avoidance for MAVs in Dynamic Environments"
% by Zhu and Alonso-Mora

function [u_opt, x_pred] = prob_collision_avoidance(x0, goal, obstacles, other_robots, params)
    % Implement probabilistic collision avoidance for 2D car model
    % 
    % Inputs:
    %   x0: Initial state [px, py, psi]
    %   goal: Goal position [px, py]
    %   obstacles: Cell array of obstacles, each with fields:
    %     - pos: Mean position [px, py]
    %     - cov: Covariance matrix (2x2)
    %     - shape: 'circle' or 'ellipse'
    %     - radius: For circle
    %     - semi_axes: [a, b] for ellipse
    %     - orientation: Orientation angle (rad) for ellipse
    %   other_robots: Cell array of other robots, each with fields:
    %     - pos: Mean position [px, py]
    %     - cov: Covariance matrix (2x2)
    %     - radius: Robot radius
    %   params: Structure with fields:
    %     - dt: Time step
    %     - N: Prediction horizon
    %     - delta_r: Collision probability threshold for robots
    %     - delta_o: Collision probability threshold for obstacles
    %     - Q: State cost matrix
    %     - R: Control cost matrix
    %     - v_min, v_max: Velocity constraints
    %     - omega_min, omega_max: Angular velocity constraints
    %     - radius: Radius of the robot
    %
    % Outputs:
    %   u_opt: Optimal control inputs [v, omega]
    %   x_pred: Predicted state trajectory
    
    % Extract parameters
    dt = params.dt;
    N = params.N;
    delta_r = params.delta_r;
    delta_o = params.delta_o;
    radius = params.radius;
    
    % Define state and control dimensions
    nx = 3; % [px, py, psi]
    nu = 2; % [v, omega]
    
    % Define optimization variables
    x = sdpvar(nx, N+1, 'full');
    u = sdpvar(nu, N, 'full');
    
    % Set initial state constraint
    x(:,1) = x0;
    
    u(:,1) = [0,0]';
    
    % Define objective function
    objective = 0;
    
     % initialize constraints
    constraints = [];
    
    % System dynamics (unicycle model)
    for k = 1:N
        % Discretized unicycle model
        x_next = car_dynamics(x(:,k), u(:,k), dt);
%         constraints = [constraints, x(:,k+1) = x_next];
        x(:,k+1) = x_next;
        
        % Input constraints
        constraints = [constraints, params.v_min <= u(1,k) <= params.v_max];
        constraints = [constraints, params.omega_min <= u(2,k) <= params.omega_max];
        
        % Cost for control input
        objective = objective + u(:,k)'*params.R*u(:,k);
        
        % Cost for state (distance to goal)
        state_cost = (x(1:2,k) - goal)'*params.Q*(x(1:2,k) - goal);
        objective = objective + state_cost;
    end
    
    % Terminal cost
    objective = objective + (x(1:2,N+1) - goal)'*(params.Q*5)*(x(1:2,N+1) - goal);
    
    % Error function inverse for probability thresholds
    erf_inv_r = erfinv(1 - 2*delta_r);
    erf_inv_o = erfinv(1 - 2*delta_o);
    
    % Uncertainty propagation matrices
    % For simplicity, we use a simplified uncertainty model
    % where the uncertainty grows with time
    Sigma_i = zeros(2, 2, N+1);
    Sigma_i(:,:,1) = [0.01, 0; 0, 0.01]; % Initial position uncertainty
    
    for k = 2:N+1
        % Simplified uncertainty propagation (grows with time)
        Sigma_i(:,:,k) = Sigma_i(:,:,k-1) + diag([0.005, 0.005]);
    end
    
    % Add collision avoidance constraints
    for k = 2:N+1
        % Robot-robot collision avoidance
        for j = 1:length(other_robots)
            robot_j = other_robots{j};
            
            % Get predicted position of robot j (assuming constant velocity for simplicity)
            % In practice, this should be from another robot's plan or prediction
            p_j = robot_j.pos;
            
            % Extract robot i's predicted position
            p_i = x(1:2,k);
            
            % Define constraint parameters as per the paper
            diff = p_i - p_j;
            dist = norm(diff);
            
            % Skip if robots are far apart
            if dist > (radius + robot_j.radius + 3*sqrt(trace(Sigma_i(:,:,k)) + trace(robot_j.cov)))
                continue;
            end
            
            % Normalized vector pointing from j to i
            a_ij = diff / max(dist, 0.01); % Avoid division by zero
            
            % Sum of radii
            b_ij = radius + robot_j.radius;
            
            % Combined covariance
            Sigma_ij = Sigma_i(:,:,k) + robot_j.cov;
            
            % Probabilistic constraint as per equation (9) in the paper
            safety_distance = erf_inv_r * sqrt(2*a_ij'*Sigma_ij*a_ij);
            constraints = [constraints, a_ij'*(p_i - p_j) - b_ij >= safety_distance];
        end
        
        % Robot-obstacle collision avoidance
        for o = 1:length(obstacles)
            obstacle = obstacles{o};
            
            % Extract robot i's predicted position
            p_i = x(1:2,k);
            p_o = obstacle.pos;
            
            if obstacle.shape == "circle"
                % For circular obstacles (simpler case)
                diff = p_i - p_o;
                dist = norm(diff);
                
                % Skip if far apart
                if dist > (radius + obstacle.radius + 3*sqrt(trace(Sigma_i(:,:,k)) + trace(obstacle.cov)))
                    continue;
                end
                
                % Normalized vector
                a_io = diff / max(dist, 0.01); % Avoid division by zero
                
                % Distance
                b_io = radius + obstacle.radius;
                
                % Combined covariance
                Sigma_io = Sigma_i(:,:,k) + obstacle.cov;
                
                % Probabilistic constraint
                safety_distance = erf_inv_o * sqrt(2*a_io'*Sigma_io*a_io);
                constraints = [constraints, a_io'*(p_i - p_o) - b_io >= safety_distance];
                
            elseif obstacle.shape == "ellipse"
                % For elliptical obstacles (more complex case as per equation 12)
                
                % Transformation matrix for the ellipse
                theta = obstacle.orientation;
                R_o = [cos(theta), -sin(theta); sin(theta), cos(theta)];
                S_o = diag(1./obstacle.semi_axes.^2);
                Omega_io = R_o'*S_o*R_o;
                Omega_sqrt = sqrtm(Omega_io); % Square root of matrix
                
                % Transform to ellipse space
                p_i_tilde = Omega_sqrt * p_i;
                p_o_tilde = Omega_sqrt * p_o;
                diff_tilde = p_i_tilde - p_o_tilde;
                dist_tilde = norm(diff_tilde);
                
                % Skip if far apart (in transformed space)
                if dist_tilde > (1 + 3*sqrt(trace(Omega_sqrt*Sigma_i(:,:,k)*Omega_sqrt')))
                    continue;
                end
                
                % Normalized vector in transformed space
                a_io_tilde = diff_tilde / max(dist_tilde, 0.01);
                
                % In transformed space, the boundary is at distance 1
                b_io = 1 + radius*norm(Omega_sqrt*[1;1])/sqrt(2); % Approximation for transformed radius
                
                % Combined covariance in transformed space
                Sigma_io_tilde = Omega_sqrt * Sigma_i(:,:,k) * Omega_sqrt';
                
                % Probabilistic constraint
                safety_distance = erf_inv_o * sqrt(2*a_io_tilde'*Sigma_io_tilde*a_io_tilde);
                constraints = [constraints, a_io_tilde'*(p_i_tilde - p_o_tilde) - b_io >= safety_distance];
            end
        end
    end
    
    % Solve optimization problem
    options = sdpsettings('verbose', 0, 'solver', 'ipopt');
    solution = optimize(constraints, objective, options);
    
    % Extract optimal control and state trajectory
    if solution.problem == 0
        u_opt = value(u);
        x_pred = value(x);
    else
        disp(['Error: ', solution.info]);
        u_opt = zeros(nu, N);
        x_pred = zeros(nx, N+1);
    end
end

function x_next = car_dynamics(x, u, dt)
    % Discrete-time unicycle model dynamics for 2D car
    % x = [px, py, psi]
    % u = [v, omega]
    
    px = x(1);
    py = x(2);
    psi = x(3);
    
    v = u(1);
    omega = u(2);
    
    % Discretized unicycle dynamics
    x_next = zeros(3, 1);
    x_next(1) = px + dt * v * cos(psi);
    x_next(2) = py + dt * v * sin(psi);
    x_next(3) = psi + dt * omega;
end

%% Example usage
function example_usage()
    % Setup parameters
    params = struct();
    params.dt = 0.1;           % Time step (s)
    params.N = 10;             % Prediction horizon
    params.delta_r = 0.03;     % Collision probability threshold for robots (3%)
    params.delta_o = 0.03;     % Collision probability threshold for obstacles (3%)
    params.Q = diag([10, 10]); % State cost matrix
    params.R = diag([1, 0.1]); % Control cost matrix
    params.v_min = -2;         % Minimum velocity (m/s)
    params.v_max = 2;          % Maximum velocity (m/s)
    params.omega_min = -1;     % Minimum angular velocity (rad/s)
    params.omega_max = 1;      % Maximum angular velocity (rad/s)
    params.radius = 0.3;       % Robot radius (m)
    
    % Initial state and goal
    x0 = [0; 0; 0];            % Initial state [px, py, psi]
    goal = [5; 5];             % Goal position [px, py]
    
    % Define obstacles
    obstacles = {};
    obstacles{1} = struct('pos', [2; 2], 'cov', diag([0.01, 0.01]), ...
                         'shape', 'circle', 'radius', 0.5);
    obstacles{2} = struct('pos', [3; 4], 'cov', diag([0.02, 0.01]), ...
                         'shape', 'ellipse', 'semi_axes', [0.6, 0.3], 'orientation', pi/4);
    
    % Define other robots
    other_robots = {};
    other_robots{1} = struct('pos', [4; 1], 'cov', diag([0.01, 0.01]), 'radius', 0.3);
    
    % Run the collision avoidance algorithm
    [u_opt, x_pred] = prob_collision_avoidance(x0, goal, obstacles, other_robots, params);
    
    % Visualize the results
    visualize_results(x0, goal, obstacles, other_robots, x_pred, params);
end

function visualize_results(x0, goal, obstacles, other_robots, x_pred, params)
    figure;
    hold on;
    
    % Plot trajectory
    plot(x_pred(1,:), x_pred(2,:), 'b-', 'LineWidth', 2);
    
    % Plot initial position and goal
    plot(x0(1), x0(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(goal(1), goal(2), 'g*', 'MarkerSize', 10);
    
    % Plot obstacles
    for i = 1:length(obstacles)
        obstacle = obstacles{i};
        if obstacle.shape == "circle"
            viscircles(obstacle.pos', obstacle.radius, 'Color', 'r');
            
            % Plot uncertainty ellipse (95% confidence)
            [V, D] = eig(obstacle.cov);
            theta = linspace(0, 2*pi, 100);
            ellipse = 2.447 * V * sqrtm(D) * [cos(theta); sin(theta)];
            plot(obstacle.pos(1) + ellipse(1,:), obstacle.pos(2) + ellipse(2,:), 'r--');
        elseif obstacle.shape == "ellipse"
            theta = linspace(0, 2*pi, 100);
            R = [cos(obstacle.orientation), -sin(obstacle.orientation); 
                 sin(obstacle.orientation), cos(obstacle.orientation)];
            ellipse = R * diag(obstacle.semi_axes) * [cos(theta); sin(theta)];
            fill(obstacle.pos(1) + ellipse(1,:), obstacle.pos(2) + ellipse(2,:), 'r', 'FaceAlpha', 0.3);
            
            % Plot uncertainty ellipse (95% confidence)
            [V, D] = eig(obstacle.cov);
            unc_ellipse = 2.447 * V * sqrtm(D) * [cos(theta); sin(theta)];
            plot(obstacle.pos(1) + unc_ellipse(1,:), obstacle.pos(2) + unc_ellipse(2,:), 'r--');
        end
    end
    
    % Plot other robots
    for i = 1:length(other_robots)
        robot = other_robots{i};
        viscircles(robot.pos', robot.radius, 'Color', 'b');
        
        % Plot uncertainty ellipse (95% confidence)
        [V, D] = eig(robot.cov);
        theta = linspace(0, 2*pi, 100);
        ellipse = 2.447 * V * sqrtm(D) * [cos(theta); sin(theta)];
        plot(robot.pos(1) + ellipse(1,:), robot.pos(2) + ellipse(2,:), 'b--');
    end
    
    % Plot robot footprint at different time steps
    robot_radius = params.radius;
    for k = 1:2:size(x_pred, 2)
        theta = linspace(0, 2*pi, 20);
        x_robot = robot_radius * cos(theta);
        y_robot = robot_radius * sin(theta);
        
        % Rotate and translate according to robot pose
        R = [cos(x_pred(3,k)), -sin(x_pred(3,k)); 
             sin(x_pred(3,k)), cos(x_pred(3,k))];
        robot_boundary = R * [x_robot; y_robot];
        plot(x_pred(1,k) + robot_boundary(1,:), x_pred(2,k) + robot_boundary(2,:), 'b-', 'LineWidth', 0.5);
        
        % Show orientation
        line([x_pred(1,k), x_pred(1,k) + robot_radius*cos(x_pred(3,k))], ...
             [x_pred(2,k), x_pred(2,k) + robot_radius*sin(x_pred(3,k))], 'Color', 'b');
    end
    
    title('Probabilistic Collision Avoidance for 2D Car Model');
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;
    grid on;
    legend('Trajectory', 'Start', 'Goal', 'Location', 'best');
end