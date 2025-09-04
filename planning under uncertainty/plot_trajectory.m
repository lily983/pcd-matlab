function plot_trajectory(d, obstacles)
    % Create a new figure
    figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Set up the plot
    axis equal;
    grid on;
    
    % Plot obstacles
    for i = 1:length(obstacles)
        obstacle = obstacles{i};
        pos = obstacle.pos;
        radius = obstacle.radius;
        
        % Plot the obstacle as a filled circle
        theta = linspace(0, 2*pi, 100);
        x_circle = pos(1) + radius * cos(theta);
        y_circle = pos(2) + radius * sin(theta);
        fill(x_circle, y_circle, 'r', 'FaceAlpha', 0.5);
        
        % Add uncertainty ellipse based on covariance matrix
        if isfield(obstacle, 'cov') && ~isempty(obstacle.cov)
            plot_covariance_ellipse(pos, obstacle.cov, 0.95, 'r--');
        end
    end
    
    % Get actual trajectory length (in case simulation ended early)
    valid_steps = find(~isnan(d.s.x(1,:)), 1, 'last');
    if isempty(valid_steps)
        valid_steps = size(d.s.x, 2);
    end
    
    % Plot the robot trajectory
    plot(d.s.x(1, 1:valid_steps), d.s.x(2, 1:valid_steps), 'b-', 'LineWidth', 2);
    
    % Mark start and goal positions
    plot(d.s.x(1, 1), d.s.x(2, 1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(d.p.goal(1), d.p.goal(2), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Plot the robot shape at regular intervals
    step_size = max(1, floor(valid_steps/15)); % Show at most 15 robot instances
    for i = 1:step_size:valid_steps
        plot_robot_shape(d.s.x(1:3, i), d.p.radius, i, valid_steps);
    end
    
    % Add legend and labels
    legend('Obstacle', 'Uncertainty', 'Robot Trajectory', 'Start', 'Goal', 'Robot');
    xlabel('X Position');
    ylabel('Y Position');
    title('Robot Trajectory with Obstacles');
    
    % Add text annotations for start and goal
    text(d.s.x(1, 1) + 0.2, d.s.x(2, 1) + 0.2, 'Start');
    text(d.p.goal(1) + 0.2, d.p.goal(2) + 0.2, 'Goal');
    
    hold off;
end

function plot_robot_shape(state, radius, step_index, total_steps)
    % Extract position and orientation
    x = state(1);
    y = state(2);
    
    % Determine if state is in radians or degrees
    % Look at the magnitude - if it's typically small (e.g., < 10), assume radians
    % Otherwise, assume degrees
    theta = state(3);
    if abs(theta) < 10
        % Likely already in radians
        theta_rad = theta;
    else
        % Likely in degrees, convert to radians
        theta_rad = theta * pi/180;
    end
    
    % Create a triangular robot shape
    robot_shape = create_robot_shape(x, y, theta_rad, radius);
    
    % Calculate color based on position in trajectory (blue → green gradient)
    color_val = step_index / total_steps;
    robot_color = [0, 0.6 + 0.4*color_val, 0.8 - 0.8*color_val];
    
    % Fill the robot shape
    fill(robot_shape.x, robot_shape.y, robot_color, 'FaceAlpha', 0.7, 'EdgeColor', [0 0 0]);
    
    % Draw a line indicating the orientation
    head_x = x + 0.8*radius * cos(theta_rad);
    head_y = y + 0.8*radius * sin(theta_rad);
    line([x, head_x], [y, head_y], 'Color', 'k', 'LineWidth', 1.5);
    
    % Add time step label for longer simulations if needed
    if total_steps > 30
        text(x + 0.2*radius, y + 0.2*radius, num2str(step_index), 'FontSize', 8);
    end
end

function robot = create_robot_shape(x, y, theta, radius)
    % Create a triangular robot with a rounded back
    
    % Front point (heading direction)
    front_x = x + radius * cos(theta);
    front_y = y + radius * sin(theta);
    
    % Back left and right points (at angle theta ± 2.5)
    back_angle_left = theta + 3*pi/4;
    back_angle_right = theta - 3*pi/4;
    
    back_left_x = x + 0.8*radius * cos(back_angle_left);
    back_left_y = y + 0.8*radius * sin(back_angle_left);
    
    back_right_x = x + 0.8*radius * cos(back_angle_right);
    back_right_y = y + 0.8*radius * sin(back_angle_right);
    
    % Create smooth back using multiple points
    n_back_points = 10;
    back_angles = linspace(back_angle_left, back_angle_right + 2*pi, n_back_points);
    back_x = zeros(1, n_back_points);
    back_y = zeros(1, n_back_points);
    
    for i = 1:n_back_points
        back_x(i) = x + 0.7*radius * cos(back_angles(i));
        back_y(i) = y + 0.7*radius * sin(back_angles(i));
    end
    
    % Combine all points to create the robot shape
    robot.x = [front_x, back_x, front_x];
    robot.y = [front_y, back_y, front_y];
end

function plot_covariance_ellipse(pos, cov, confidence, line_style)
    % Calculate eigenvalues and eigenvectors of covariance matrix
    [eigvec, eigval] = eig(cov);
    
    % Calculate the scaling factor for the desired confidence level
    % For 2D Gaussian: 95% confidence = 5.991 chi-squared value
    if confidence == 0.95
        s = 5.991;
    elseif confidence == 0.99
        s = 9.21;
    else
        s = 2.447; % 68% confidence
    end
    
    % Scale eigenvalues for desired confidence ellipse
    eigval = sqrt(s * eigval);
    
    % Generate points for ellipse
    theta = linspace(0, 2*pi, 100);
    ellipse = zeros(2, length(theta));
    for i = 1:length(theta)
        ellipse(:, i) = pos + eigvec * [eigval(1,1)*cos(theta(i)); eigval(2,2)*sin(theta(i))];
    end
    
    % Plot the ellipse
    plot(ellipse(1, :), ellipse(2, :), line_style, 'LineWidth', 1.5);
end