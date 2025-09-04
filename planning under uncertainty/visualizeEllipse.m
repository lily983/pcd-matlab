function visualizeEllipse(a, b, tc, theta)
% visualizeEllipse plots an ellipse given its semi-axes, center, and rotation angle.
%
%   Inputs:
%       a     - Semi-axis length along the x-direction (before rotation)
%       b     - Semi-axis length along the y-direction (before rotation)
%       tc    - 1x2 vector specifying the ellipse center [tc_x, tc_y]
%       theta - Rotation angle (in radians)
%
%   Example:
%       visualizeEllipse(5, 3, [2, -1], pi/4);

    % Parameterize the ellipse using the angle t from 0 to 2*pi
    t = linspace(0, 2*pi, 100);
    
    % Parametric equations of the rotated ellipse
    % These equations come from:
    % x = tc_x + a*cos(t)*cos(theta) - b*sin(t)*sin(theta)
    % y = tc_y + a*cos(t)*sin(theta) + b*sin(t)*cos(theta)
    x = tc(1) + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
    y = tc(2) + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);
    
    % Create the plot
%     figure;
    plot(x, y, 'color', hex2rgb('45AC59'), 'LineWidth', 2);
%     hold on;
%     plot(tc(1), tc(2), 'ro', 'MarkerFaceColor', 'r'); % plot the center of the ellipse
%     hold off;
    
    % Set plot properties
%     axis equal;  % Ensure equal scaling for x and y axes
%     grid on;
%     xlabel('X');
%     ylabel('Y');
%     title('Rotated Ellipse Visualization');
end
