function plot_trajectory_with_uncertainty(d, obstacles)
    % Create a new figure
    figure('Position', [100, 100, 900, 700]);
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
    % Find the last non-zero entry - if there are zeros at the end of the trajectory
    nonzero_entries = find(abs(d.s.x(1,:)) > 1e-10 | abs(d.s.x(2,:)) > 1e-10);
    if ~isempty(nonzero_entries)
        valid_steps = nonzero_entries(end);
    else
        valid_steps = 1; % At least include the starting position
    end
    
    % Plot the robot trajectory
    plot(d.s.x(1, 1:valid_steps), d.s.x(2, 1:valid_steps), 'b-', 'LineWidth', 2);
    
    % Mark start and goal positions
    plot(d.s.x(1, 1), d.s.x(2, 1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(d.p.goal(1), d.p.goal(2), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Initialize uncertainty covariance with initial value
    Sigma = d.p.Sigma_0;
    
    % Plot the robot and its uncertainty at regular intervals
    step_size = max(1, floor(valid_steps/10)); % Show at most 10 robot instances with uncertainty
    for i = 1:step_size:valid_steps
        % Plot robot
%         plot_circular_robot(d.s.x(1:3, i), d.p.radius, i, valid_steps);
        
        % Plot position uncertainty as ellipse
        pos = d.s.x(1:2, i);
        
        % For position uncertainty, use only the position part of the covariance
        % Assuming Sigma contains position covariance in the first 2x2 block
        pos_cov = Sigma(1:2, 1:2);
        
        % Plot the uncertainty ellipse
        % Use a different color based on time step (from light blue to dark blue)
        color_intensity = 0.8 - 0.6 * (i / valid_steps);
        ellipse_color = [color_intensity, color_intensity, 0.8];
        
        plot_colored_covariance_ellipse(pos, pos_cov, 0.95, ellipse_color);
        
        % Propagate uncertainty for the next step
        % Sigma_i = Sigma_i + d.p.Sigma_propogagtion
        if i < valid_steps
            Sigma = Sigma + d.p.Sigma_propogagtion;
        end
    end
    
    % Add legend and labels
    legend('Obstacle', 'Obstacle Uncertainty', 'Robot Trajectory', 'Start', 'Goal', 'Robot', 'Position Uncertainty');
    xlabel('X Position');
    ylabel('Y Position');
    title('Robot Trajectory with Propagated Position Uncertainty');
    
    % Add text annotations for start and goal
    text(d.s.x(1, 1) + 0.2, d.s.x(2, 1) + 0.2, 'Start');
    text(d.p.goal(1) + 0.2, d.p.goal(2) + 0.2, 'Goal');
    
    % Add colorbar to indicate uncertainty growth over time
    colormap(flipud(winter));
    c = colorbar;
    c.Label.String = 'Time Step (Uncertainty Growth)';
    c.Ticks = linspace(0, 1, 5);
    c.TickLabels = {'Later Steps', '', 'Middle Steps', '', 'Early Steps'};
    
    hold off;
end

function plot_circular_robot(state, radius, step_index, total_steps)
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
    
    % Calculate color based on position in trajectory (blue → green gradient)
    color_val = step_index / total_steps;
    robot_color = [0, 0.6 + 0.4*color_val, 0.8 - 0.8*color_val];
    
    % Plot the robot as a circle
    theta_circle = linspace(0, 2*pi, 50);
    x_circle = x + radius * cos(theta_circle);
    y_circle = y + radius * sin(theta_circle);
    fill(x_circle, y_circle, robot_color, 'FaceAlpha', 0.5, 'EdgeColor', [0 0 0]);
    
    % Draw a line indicating the orientation
    head_x = x + radius * cos(theta_rad);
    head_y = y + radius * sin(theta_rad);
    line([x, head_x], [y, head_y], 'Color', 'k', 'LineWidth', 1.5);
    
    % Add a small dot at the center for clarity
    plot(x, y, 'k.', 'MarkerSize', 5);
    
    % Add time step label if needed for longer simulations
    if total_steps > 30 && mod(step_index, 5) == 0
        text(x + 1.2*radius, y, num2str(step_index), 'FontSize', 8);
    end
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

function plot_colored_covariance_ellipse(pos, cov, confidence, color)
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
    
    % Plot the ellipse with specified color
    plot(ellipse(1, :), ellipse(2, :), 'Color', color, 'LineWidth', 1.5);
end

% Example usage:
% d = load('robot_data.mat'); % Load your data
% obstacles = {...}; % Define obstacles
% plot_trajectory_with_uncertainty(d, obstacles);