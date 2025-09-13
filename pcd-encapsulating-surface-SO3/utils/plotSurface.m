
function plotSurface(X, N, color, alpha)
% figure; hold on; axis equal

if nargin==3
    alpha=0.3;
end
x = reshape(X(1,:), N(1), N(2));
y = reshape(X(2,:),  N(1), N(2));
z = reshape(X(3,:),  N(1), N(2));

if ~all(color >= 0 & color <= 1)
    color = color/255;
end

surf(x, y, z, 'FaceColor', color, 'EdgeColor', color,...
                'FaceAlpha', alpha, 'EdgeAlpha', alpha);
end