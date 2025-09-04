clc; clear; close all;

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

%% Loading rotations from samples
T = csv2group(summaryposedata, 'PCG3');
[T_filter, filtered_indices] = filterData(T, eye(3));

% T_filter = T;

R = T_filter(1:3,1:3,:);
tc = T_filter(1:3,4,:);
mu = get_mean_cov(R, 'SO', 0);

[mu_t, Sigma_t] = get_mean_cov(T_filter, 'R', 0);

[mu_SO3, Sigma_SO3, ~] = get_mean_cov(R, 'SO', 0);

%%
s1_tc_list = zeros(3, 100);
s2_tc_list = zeros(3, 100);

for i = 1:100
    s1_tc_list(:, i) = mean(tc,3) - 0.1*rand(3,1);
    s2_tc_list(:, i) = mean(tc,3) +  0.2*rand(3,1);
end
    
%%
prob_list = zeros(1,100);

prob_sq_list = zeros(1,100);

for i = 1:100
    s1 = SuperQuadrics({[0.2927352017397958,0.29158779787962896,0.12707043772973667]', [1.6698    0.5169], [0, 0]...
                 s1_tc_list(:, i),  rotm2quat(mu_SO3), [20, 20]});

    s2 = SuperQuadrics({[0.26966692120426106,0.18009964067697282,0.10951677121080064]',  [1.2826    1.0934], [0, 0]...
        s2_tc_list(:, i), rotm2quat(mu_SO3), [20, 20]});

    scale = 1.0;

    sub1 = EnlargedSuperQuadrics(s1, mu_SO3, scale);

    sub2 = EnlargedSuperQuadrics(s2, mu_SO3, scale);

    minkSum_sub = EnlargedMinkSumClosedForm(sub1, sub2);

    xx = s2.tc - s1.tc;
    Sigmax = Sigma_t;


    [prob, a, x_mink, a_T, x_minkT] = linearChanceConstraintEnlargedSQ(sub1, sub2, xx, Sigmax, 0);
    
    prob_sq=  linearChanceConstraintSQ(s1, s2, xx, Sigmax, 'tangent-point-cfc', 0);
    
    prob_list(i) = prob;
    
    prob_sq_list(i) = prob_sq;
    
end

%%
s1 = SuperQuadrics({[0.2927352017397958,0.29158779787962896,0.12707043772973667]', [1.6698    0.5169], [0, 0]...
             mean(tc,3) - [0.1;0.06;0.09],  rotm2quat(mu_SO3), [20, 20]});
        
s2 = SuperQuadrics({[0.26966692120426106,0.18009964067697282,0.10951677121080064]',  [1.2826    1.0934], [0, 0]...
    mean(tc,3) + [0.2;0.2;0.2], rotm2quat(mu_SO3), [20, 20]});

scale = 1.2;

sub1 = EnlargedSuperQuadrics(s1, R, scale);

sub2 = EnlargedSuperQuadrics(s2, R, scale);

minkSum_sub = EnlargedMinkSumClosedForm(sub1, sub2);

xx = s2.tc - s1.tc;
Sigmax = Sigma_t;

%%
[prob, a, x_mink, a_T, x_minkT] = linearChanceConstraintEnlargedSQ(sub1, sub2, xx, Sigmax, 0);

% figure; hold on; visualizePositionErrors(s1, s2,  sigma,  sigma, s1_color, s2_color)

%% Plot figures

% Visualize superquadrics
figure; hold on; axis equal
s1.PlotShape(s1_color, 0.3);
for i=1:size(R,3)
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    s1_i = SuperQuadrics({s1.a, s1.eps, [0,0], tc(1:3,i) - [0.1;0.06;0.09], rotm2quat(R(:,:,i)), [20,20]});
    s1_i.PlotShape(s1_color, 0.2,0.2);
    pause(0.1);
end
sub1.PlotShape(s2_color, 0.8);
% 
figure; hold on; axis equal
s2.PlotShape(s2_color, 0.3);
for i=1:size(R,3)
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    s2_i = SuperQuadrics({s2.a, s2.eps, [0,0], tc(1:3,i) + [0.2;0.2;0.2], rotm2quat(R(:,:,i)), [20,20]});
    s2_i.PlotShape(s2_color, 0.2,0.2);
    pause(0.1);
end
sub2.PlotShape(s1_color, 0.8);

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

%% Plot figures

% Visualize superquadrics
figure; hold on; axis equal

sub1.PlotShape(s1_color, 0.3);
for i=1:size(R,3)
    s1_i = SuperQuadrics({s1.a, s1.eps, [0,0], tc(1:3,i) - [0.1;0.06;0.09], [1 0 0 0], [20,20]});
    sub1_i = EnlargedSuperQuadrics(s1_i, R, scale);
    sub1_i.PlotShape(s1_color, 0.2,0.2);
    pause(0.1);
end
% sub1.PlotShape(s1_color, 0.3);
% 
% figure; hold on; axis equal
sub2.PlotShape(s2_color, 0.3);
for i=1:size(R,3)
    s2_i = SuperQuadrics({s2.a, s2.eps, [0,0], tc(1:3,i) + [0.2;0.2;0.2], [1 0 0 0], [20,20]});
    sub2_i = EnlargedSuperQuadrics(s2_i, R, scale);
    sub2_i.PlotShape(s2_color, 0.2,0.2);
    pause(0.1);
end

xLimits = [-0.4,0.8];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.5

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 0.05; % Control the size of the arrowhead

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

% print('sub1sub2positionerror', '-dpng', '-r600')

%% Visualize the Minkowski sum and confidence ellipsoid of the position error
figure; hold on

minksum_sub = EnlargedMinkSumClosedForm(sub1, sub2);

minksum_sub.PlotShape(mink_sf_color,0.5,1)

% plot position error ellipse
scatter3(xx(1), xx(2), xx(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 5);
% plot_ellipse(Sigmax, xx, 'confidence', 0.5, 'edgecolor', position_error_color)
plot_ellipse(Sigmax, xx, 'confidence', 7, ...
    'edgecolor', position_error_color, 'fillcolor', position_error_color, 'alpha', 0.3)

axis off;
axis equal
xLimits = [-0.2,0.6];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.5

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 0.05; % Control the size of the arrowhead

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


%% Visualize the Minkowski sum and confidence ellipsoid of the position error
figure; hold on

% Space transformation so that Sigma_x = eye(3)
xx_T = Sigmax^0.5 \ xx;

x_mink = minksum_sub.GetPoints();
x_mink_T = Sigmax^0.5 \ x_mink;
N = s1.N;
X_T = reshape(x_mink_T(1,:), N(1), N(2)); %change array to matrix form
Y_T = reshape(x_mink_T(2,:), N(1), N(2));
Z_T = reshape(x_mink_T(3,:), N(1), N(2));

mink_sf_T = surf(X_T, Y_T, Z_T,...
 'EdgeColor',mink_sf_color, 'EdgeAlpha', 1,...
 'FaceAlpha', 0.2, 'FaceColor', mink_sf_color); %plot contact space mud yellow


% plot position error ellipse
scatter3(xx_T(1), xx_T(2), xx_T(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);
% plot_ellipse(Sigmax, xx, 'confidence', 0.5, 'edgecolor', position_error_color)
% plot_ellipse(eye(3), xx_T, 'confidence', 5, ...
%     'edgecolor', position_error_color, 'fillcolor', position_error_color, 'alpha', 0.3)

[prob, a, x_mink, a_T, x_minkT] = linearChanceConstraintEnlargedSQ(sub1, sub2, xx, Sigmax, 0);
% plotPlane(a_T, x_minT, tangent_plane_color);

scatter3(x_minkT(1), x_minkT(2), x_minkT(3), 'MarkerFaceColor', tangent_plane_color,...
 'MarkeredgeColor', tangent_plane_color, 'SizeData', 20);

plot3([x_minkT(1), xx_T(1)], [x_minkT(2), xx_T(2)], [x_minkT(3), xx_T(3)], '-', 'LineWidth', 2, 'Color', tangent_plane_color);

axis off;
axis equal
xLimits = [-20, 20];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.1

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 2; % Control the size of the arrowhead

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

% print('transformedminkpositionerror', '-dpng', '-r600')


%% Visualize the Minkowski sum and confidence ellipsoid of the position error and tangent plane
figure; hold on

% minksum_sub = EnlargedMinkSumClosedForm(sub1, sub2);

minksum_sub.PlotShape(mink_sf_color,0.5,1)

% plot position error ellipse
scatter3(xx(1), xx(2), xx(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);
% plot_ellipse(Sigmax, xx, 'confidence', 0.5, 'edgecolor', position_error_color)
plot_ellipse(Sigmax, xx, 'confidence', 8.2, ...
    'edgecolor', position_error_color, 'fillcolor', position_error_color, 'alpha', 0.2)

plotPlane(a, x_mink, tangent_plane_color);

axis off;
axis equal
xLimits = [-0.4,1];
yLimits = xLimits;
zLimits = xLimits;
Linewidth = 0.1

% Plot the X, Y, and Z axes at the origin with arrowheads
plot3(xLimits, [0 0], [0 0], 'k', 'LineWidth', Linewidth); % X-axis
plot3([0 0], yLimits, [0 0], 'k', 'LineWidth', Linewidth); % Y-axis
plot3([0 0], [0 0], zLimits, 'k', 'LineWidth', Linewidth); % Z-axis

% Create arrowheads using patch
arrowSize = 0.05; % Control the size of the arrowhead

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
