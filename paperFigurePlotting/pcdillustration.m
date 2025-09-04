clc; clear; close all;

% %% Generate randome superquadrics, their Minkowski sum, and bounding ellipsoids
% s1 = SuperQuadrics({1.5 * rand(1, 3) + 0.5, 2 * rand(1, 2) + 0.004, [0, 0]...
%             2 *rand(3,1), getRandomQuaternion(), [20, 20]});
%         
% s2 = SuperQuadrics({1.5 * rand(1, 3) + 0.5,  2 * rand(1, 2) + 0.004, [0, 0]...
%     2 * rand(3,1) + 2, getRandomQuaternion(),[20, 20]});
%%
s1 = SuperQuadrics({[1.4930    1.1242    1.7629], [1.6698    0.5169], [0, 0]...
            [1.2269    1.1645    1.0815]', [0.8962    0.2728    0.3277    0.1228], [20, 20]});
        
s2 = SuperQuadrics({[1.9097    1.4683    1.2192],  [1.2826    1.0934], [0, 0]...
    [3.2946    3.0878    3.4421]', [0.4549    0.8651    0.1904    0.0921],[20, 20]});

% Minkowski sum 
% Plot Minkowski Sum
m1 = s1.GetGradientsCanonical();
minkSum = MinkSumClosedForm(s1, s2, quat2rotm(s1.q), quat2rotm(s2.q));
x_mink = minkSum.GetMinkSumFromGradient(m1);
N = s1.N;
X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

%%
% Position error distribution
sigma = [0.10, 0, 0; 0, 0.10, 0; 0, 0, 0.35];
cov1 =    quat2rotm(s1.q) * sigma * quat2rotm(s1.q)';
cov2 =    quat2rotm(s2.q) * sigma * quat2rotm(s2.q)';

% %% Get the lcc-tangent
xx = s2.tc - s1.tc;
Sigmax = cov1 + cov2;
 [p, ~, a, x_mink_opt] = linearChanceConstraintSQ(s1, s2, xx, Sigmax, 'tangent-point-cfc', 0);
 
%% Get bounding ellipsoid for s1 and s2
s1Points = s1.GetPoints();
s2Points = s2.GetPoints();
[invSigma1, e1] = MinVolEllipse(s1Points, 0.001);
[invSigma2, e2] = MinVolEllipse(s2Points, 0.001);
Sigma1 = inv(invSigma1);
Sigma2 = inv(invSigma2);

%% Color 
%     '66c2a5' 'mint green'
%     '93a5cd' 'violet blue'
%     'a0cbe1' 'sky blue'
%     'e5c494' 'mud yellow'
%     'a6d854' 'chile green'
%     'fc8d62' 'orange'
%     'fed930' 'lemon'
%     'e78ac3' 'barbee pink'
%     'c8e9a0' baby green
%     '6699cc' baby blue
%     'e9b872' baby orange

s1_color = hex2rgb('66c2a5'); %mint green\
s2_color = hex2rgb('93a5cd'); %violet blue
mink_sf_color = hex2rgb('c8e9a0'); %baby green
position_error_color  = hex2rgb('e9b872'); %baby orange
tangent_plane_color = hex2rgb('6699cc'); %baby blue

%% Plot figures

% Visualize superquadrics
figure; hold on

s1.PlotShape(s1_color, 0.3);

s2.PlotShape(s2_color, 0.3);

ax = gca; % Get current axes
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis off;
axis equal

view(ax,[-17.9014170161219 19.7715123087351]);
hold(ax,'off');
% Set the remaining axes properties
set(ax,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
    [1.61040886050011 1.27014505287831 1],'XAxisLocation','origin',...
    'YAxisLocation','origin');


%% Visualize the effect of the position errors
figure; hold on
s1.PlotShape(s1_color, 0.2, 0.8);
s2.PlotShape(s2_color, 0.2, 0.8);

visualizePositionErrors(s1, s2,  sigma,  sigma, s1_color, s2_color)

axis off;
axis equal
xLimits = [-2,2];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.1

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 0.2; % Control the size of the arrowhead

% X-axis arrowhead
patch([xLimits(2) xLimits(2)-arrowSize xLimits(2)-arrowSize], ...
      [0 -arrowSize/2 arrowSize/2], ...
      [0 0 0], 'k');

% Y-axis arrowhead
patch([0 0 0], ...
      [yLimits(2) yLimits(2)-arrowSize yLimits(2)-arrowSize], ...
      [0 -arrowSize/2 arrowSize/2], 'k');

% Z-axis arrowhead
patch([-arrowSize/2 arrowSize/2 0], ...
      [0 0 0], ...
      [zLimits(2)-arrowSize zLimits(2)-arrowSize zLimits(2)], 'k');

axes1 = gca; 
view(axes1,[121.854311694003 21.6347595905472]);
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
    [1.2678936605317 1 1.18777820973877]);
% Create textarrow
annotation('textarrow',[0.68,0.65],[0.49,0.56],...
    'String','position error $\mathbf{x}_2$',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.56,0.54],[0.38,0.43],...
    'String','position error $\mathbf{x}_1$',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.45 0.48],[0.65 0.59],'String','$S_1^{UB}$',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.5 0.53],[0.65 0.57],'String','$S_2^{UB}$',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');
% print('s1s2positionerror', '-dpng', '-r600')

%% Visualize the Minkowski sum and confidence ellipsoid of the position error
figure; hold on

mink_sf = surf(X, Y, Z,...
 'EdgeColor', mink_sf_color, 'EdgeAlpha', 1,...
 'FaceAlpha', 0.5, 'FaceColor', mink_sf_color); %plot contact space mud yellow

% plot position error ellipse
scatter3(xx(1), xx(2), xx(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 5);
% plot_ellipse(Sigmax, xx, 'confidence', 0.5, 'edgecolor', position_error_color)
plot_ellipse(Sigmax, xx, 'confidence', 1, ...
    'edgecolor', position_error_color, 'fillcolor', position_error_color, 'alpha', 0.3)

axis off;
axis equal
xLimits = [-2,2];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.1

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 0.2; % Control the size of the arrowhead

% X-axis arrowhead
patch([xLimits(2) xLimits(2)-arrowSize xLimits(2)-arrowSize], ...
      [0 -arrowSize/2 arrowSize/2], ...
      [0 0 0], 'k');

% Y-axis arrowhead
patch([0 0 0], ...
      [yLimits(2) yLimits(2)-arrowSize yLimits(2)-arrowSize], ...
      [0 -arrowSize/2 arrowSize/2], 'k');

% Z-axis arrowhead
patch([-arrowSize/2 arrowSize/2 0], ...
      [0 0 0], ...
      [zLimits(2)-arrowSize zLimits(2)-arrowSize zLimits(2)], 'k');

axes1 = gca; 
view(axes1,[121.854311694003 21.6347595905472]);
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
    [1.2678936605317 1 1.18777820973877]);
% Create textarrow
annotation('textarrow',[0.7,0.62],[0.32,0.39],...
    'String',{'Minkowski sum boundary', '$\mathbf{x}_\Sigma\in S_1^{UB}\oplus -S_2^{UB}$'},...
    'HorizontalAlignment','center',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.63,0.64],[0.77,0.73],...
    'String',{'relative position error','$\mathbf{x}\sim \mathcal{N}(\mathbf{p}_x, \Sigma_x)$'},...
    'HorizontalAlignment','center',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% print('minkpositionerror', '-dpng', '-r600')

%% Visualize the transformed Minkowski sum and confidence ellipsoid of the position error and tangent plane
figure; hold on;

% Space transformation so that Sigma_x = eye(3)
xx_T = Sigmax^0.5 \ xx;
minkSum_T = MinkSumClosedForm(s1, s2, Sigmax^0.5 \ quat2rotm(s1.q), Sigmax^0.5 \ quat2rotm(s2.q));
x_mink_T = Sigmax^0.5 \ x_mink;
X_T = reshape(x_mink_T(1,:), N(1), N(2)); %change array to matrix form
Y_T = reshape(x_mink_T(2,:), N(1), N(2));
Z_T = reshape(x_mink_T(3,:), N(1), N(2));

mink_sf_T = surf(X_T, Y_T, Z_T,...
 'EdgeColor',mink_sf_color, 'EdgeAlpha', 1,...
 'FaceAlpha', 0.2, 'FaceColor', mink_sf_color); %plot contact space mud yellow

 [~, ~, ~, ~, m_T, x_mink_T] = linearChanceConstraintSQ(s1, s2, xx, Sigmax, 'tangent-point-cfc', 0);
 a_T = m_T./norm(m_T);

plotPlane(a_T, x_mink_T, tangent_plane_color);

% plot position error ellipse
scatter3(xx_T(1), xx_T(2), xx_T(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 10);

plot_ellipse(eye(3)* (a_T' * xx_T - a_T' * x_mink_T)^2, xx_T, 'confidence', 1, ...
    'edgecolor', position_error_color, 'fillcolor', position_error_color, 'alpha', 0.2)

axis off;
axis equal
xLimits = [-2,2];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.1

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 0.2; % Control the size of the arrowhead

% X-axis arrowhead
patch([xLimits(2) xLimits(2)-arrowSize xLimits(2)-arrowSize], ...
      [0 -arrowSize/2 arrowSize/2], ...
      [0 0 0], 'k');

% Y-axis arrowhead
patch([0 0 0], ...
      [yLimits(2) yLimits(2)-arrowSize yLimits(2)-arrowSize], ...
      [0 -arrowSize/2 arrowSize/2], 'k');

% Z-axis arrowhead
patch([-arrowSize/2 arrowSize/2 0], ...
      [0 0 0], ...
      [zLimits(2)-arrowSize zLimits(2)-arrowSize zLimits(2)], 'k');

axes1 = gca; 
view(axes1,[121.854311694003 21.6347595905472]);
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
    [1.2678936605317 1 1.18777820973877]);
% Create textarrow
annotation('textarrow',[0.65,0.6],[0.32,0.47],...
    'String',{'transformed position error','$\mathcal{N}(\mathbf{m}_x^\prime, \textbf{I}_3)$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.43,0.48],[0.76,0.73],...
    'String',{'tangent plane','$\mathbf{ax}-b$'},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.42,0.45],[0.35,0.39],...
    'String','$\mathbf{x}_\Sigma''$',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

% Create textarrow
annotation('textarrow',[0.4,0.546],[0.44,0.482],...
    'String','$\mathbf{x}_\Sigma''(\mathbf{n_{opt}})$',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',7,...
    'FontSize',16,...
    'FontName','Times');

print('tangentplane', '-dpng', '-r600')
%% Visualize the Minkowski sum and tangent plane
figure; hold on

mink_sf = surf(X, Y, Z,...
 'EdgeColor', mink_sf_color, 'EdgeAlpha', 1,...
 'FaceAlpha', 0.5, 'FaceColor', mink_sf_color); %plot contact space

% plot position error ellipse
scatter3(xx(1), xx(2), xx(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);
% plot_ellipse(Sigmax, xx, 'confidence', 0.5, 'edgecolor', position_error_color)
plot_ellipse(Sigmax, xx, 'confidence', 1.1, ...
    'edgecolor', position_error_color, 'fillcolor', position_error_color, 'alpha', 0.25)

plotPlane(a, x_mink_opt, tangent_plane_color);

ax = gca; % Get current axes
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis on;
axis equal

view(ax,[-17.9014170161219 19.7715123087351]);
hold(ax,'off');
% Set the remaining axes properties
set(ax,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
    [1.61040886050011 1.27014505287831 1],'XAxisLocation','origin',...
    'YAxisLocation','origin');
