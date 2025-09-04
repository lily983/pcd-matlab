function visualize_results(d, obstacles)
    % Create a figure
    figure('Position', [100, 100, 1000, 800]);
    
    % Create subplot layout
    subplot(2, 2, [1, 3]); % Main trajectory plot
    hold on;
    grid on;
    
    % Plot trajectory
    trajectory = d.s.x(1:2, 1:end);
    plot(trajectory(1, :), trajectory(2, :), 'b-', 'LineWidth', 2);
    plot(trajectory(1, 1), trajectory(2, 1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Start
    plot(d.p.goal(1), d.p.goal(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Goal
    
    % Plot robot at several time steps
    robot_radius = d.p.radius;
    step_size = max(1, floor(size(d.s.x, 2) / 10)); % Plot up to 10 robot positions
    
    for i = 1:step_size:size(d.s.x, 2)
        px = d.s.x(1, i);
        py = d.s.x(2, i);
        theta = d.s.x(3, i);
        
        % Draw robot as a circle
        draw_robot(px, py, robot_radius, theta, [0.7, 0.7, 1.0]);
        
        % Add time marker
        text(px+0.1, py+0.1, sprintf('t=%0.1f', (i-1)*d.p.dt), 'FontSize', 8);
    end
    
    % Plot obstacles
    for i = 1:length(obstacles)
        obstacle = obstacles{i};
        pos = obstacle.pos;
        semi_axes = obstacle.semi_axes;
        orientation = obstacle.orientation;
        
        % Draw elliptical obstacle
        draw_ellipse(pos(1), pos(2), semi_axes(1), semi_axes(2), orientation, [1.0, 0.7, 0.7]);
    end
    
    % Add labels and title
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Robot Trajectory and Obstacles');
    axis equal;
    
    % Add grid of arrows showing avoidance field (optional)
    [X, Y] = meshgrid(linspace(min(trajectory(1,:))-0.5, max(trajectory(1,:))+0.5, 15), ...
                      linspace(min(trajectory(2,:))-0.5, max(trajectory(2,:))+0.5, 15));
    U = zeros(size(X));
    V = zeros(size(Y));
    
    % Calculate repulsive vectors
    for i = 1:numel(X)
        pos = [X(i); Y(i)];
        repulse = [0; 0];
        
        % Calculate repulsion from all obstacles
        for j = 1:length(obstacles)
            obstacle = obstacles{j};
            obs_pos = obstacle.pos;
            
            % Vector from obstacle to point
            diff = pos - obs_pos;
            dist = norm(diff) - max(obstacle.semi_axes);
            
            % Apply repulsion if within influence radius
            if dist < d.p.repulsion_radius && dist > 0
                repulse = repulse + diff/norm(diff) * d.p.obstacle_repulsion * (1 - dist/d.p.repulsion_radius)^2;
            end
        end
        
        % Add attraction to goal
        goal_diff = d.p.goal(1:2) - pos;
        goal_dist = norm(goal_diff);
        if goal_dist > 0
            goal_attraction = goal_diff/goal_dist * 0.5;
            repulse = repulse + goal_attraction;
        end
        
        % Normalize and store
        if norm(repulse) > 0
            repulse = repulse / norm(repulse) * 0.15;
        end
        
        U(i) = repulse(1);
        V(i) = repulse(2);
    end
    
    % Plot vector field
    quiver(X, Y, U, V, 0.5, 'k');
    
    % Create subplot for control inputs
    subplot(2, 2, 2);
    t = (0:size(d.s.u, 2)-1) * d.p.dt;
    plot(t, d.s.u(1, :), 'b-', 'LineWidth', 2);
    hold on;
    plot(t, d.s.u(2, :), 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Control Input');
    title('Control Inputs');
    legend('v (m/s)', '\omega (rad/s)');
    
    % Create subplot for computation time
    subplot(2, 2, 4);
    t_cpu = (1:length(d.s.CPU_time)) * d.p.dt;
    plot(t_cpu, d.s.CPU_time, 'k-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('CPU Time (s)');
    title('NMPC Computation Time');
    
    % Add statistics
    avg_cpu = mean(d.s.CPU_time);
    max_cpu = max(d.s.CPU_time);
    text(0.1, 0.9*max_cpu, sprintf('Avg CPU: %.3f s', avg_cpu), 'FontSize', 10);
    text(0.1, 0.8*max_cpu, sprintf('Max CPU: %.3f s', max_cpu), 'FontSize', 10);
    
    sgtitle('NMPC Car Control with Obstacle Avoidance', 'FontSize', 14);
end

function draw_robot(x, y, radius, theta, color)
    % Draw robot as a circle with orientation line
    t = linspace(0, 2*pi, 50);
    circle_x = x + radius * cos(t);
    circle_y = y + radius * sin(t);
    
    fill(circle_x, circle_y, color, 'EdgeColor', 'b');
    
    % Draw orientation line
    line_length = radius * 0.8;
    line_x = [x, x + line_length * cos(theta)];
    line_y = [y, y + line_length * sin(theta)];
    
    plot(line_x, line_y, 'k-', 'LineWidth', 2);
end

function draw_ellipse(x, y, a, b, theta, color)
    % Draw ellipse with semi-major axis a, semi-minor axis b, and orientation theta
    t = linspace(0, 2*pi, 100);
    
    % Ellipse in local coordinates
    X = a * cos(t);
    Y = b * sin(t);
    
    % Rotation matrix
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % Transform to global coordinates
    points = R * [X; Y];
    
    % Translate to center
    ellipse_x = points(1, :) + x;
    ellipse_y = points(2, :) + y;
    
    % Plot filled ellipse
    fill(ellipse_x, ellipse_y, color, 'EdgeColor', 'r');
end