% This file visualize lcc related concepts in 2D
clc; clear; close all; 

%% lcc-center: SQ case
figure; hold on; 

% Construct two SQ
s1 = SuperEllipse([2, 6, 0.2, 0, -19, -9, -0.5, 50]);
s2 = SuperEllipse([5, 3, 0.2, 0, -14, 5, -0.0, 50]);

s1_color = hex2rgb('45AC59'); %green
s1.PlotShape(s1_color, 1, 1);

s2_color = hex2rgb('4FAAD1'); %blue
s2.PlotShape(s2_color, 1, 1);

% Get their exact Minkowski sum boundary by the definition
% Move them to the origin before operation
s1Points = s1.GetPoints()' - s1.tc';
s2Points = -1*s2.GetPoints()' + s2.tc';
pgon1 = polyshape(s1Points(:,1), s1Points(:,2));
pgon2= polyshape(s2Points(:,1), s2Points(:,2));
mSum = minkowskiSum(pgon1, pgon2);

patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% First get an bounding ellipse for the Minkowski sum
[invSigmaf, mf] = MinVolEllipse(mSum.Vertices', 0.001);
Sigmaf = inv(invSigmaf);

% plot bounding ellipse
bounding_ellip_color = hex2rgb('f44336');
plot_ellipse(Sigmaf, mf, 'edgecolor', bounding_ellip_color);

% plot position error
Sigma1 = angle2rotm(s1.ang) * [4, 0; 0, 9]*1e-01 * angle2rotm(s1.ang)';
Sigma2 = angle2rotm(s2.ang) * [8, 0; 0, 15] * 1e-01 * angle2rotm(s2.ang)';

n = 20;
x1_rand = mvnrnd(zeros(2,1), Sigma1, n);
s1_shift = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
    s1.tc(1), s1.tc(2), s1.ang, s1.N]);
for i = 1:size(x1_rand, 1)
    s1_shift.tc = s1.tc + x1_rand(i,:)';
    s1_shift.PlotShape(s1_color, 0.15, 0.0);
end

x2_rand = mvnrnd(zeros(2,1), Sigma2, n);
s2_shift = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
    s2.tc(1), s2.tc(2), s2.ang, s2.N]);
for i = 1:size(x1_rand, 1)
    s2_shift.tc = s2.tc + x2_rand(i,:)';
    s2_shift.PlotShape(s2_color, 0.15, 0.0);
end

% Plot the relative position error
xx = s2.tc - s1.tc;
Sigmax =Sigma1+Sigma2;
xx_color = hex2rgb('0b5394');

% Plot relative position error
scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', xx_color, 'SizeData', 20);
plot_ellipse(Sigmax, xx, 'edgecolor', xx_color, 'confidence', 0.5);
plot_ellipse(Sigmax, xx, 'edgecolor', xx_color, 'confidence', 0.99);

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
axis equal

%% lcc-center: Zhu's method
% space transformation inv(Sigmaf)^0.5
figure; hold on;
        
% transformed bounding ellip
Sigmaf_T = Sigmaf^0.5 \ Sigmaf / Sigmaf^0.5;
[R_Sigmaf_T, Lamda_Sigmaf_T] = svd(Sigmaf_T);

% plot transformed bounding ellip(should be a sphere)
plot_ellipse(Sigmaf_T, 'edgecolor', bounding_ellip_color);

% transformed Minkowski sum
mSum_T=  Sigmaf^0.5 \ [mSum.Vertices(:,1), mSum.Vertices(:,2)]';
patch(mSum_T(1, :), mSum_T(2, :), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% transformed relative position error
xx_T = Sigmaf^0.5 \ xx;
Sigmax_T = Sigmaf^0.5 \ Sigmax / Sigmaf^0.5;

% plot transformed position error
scatter(xx_T(1), xx_T(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', xx_color, 'SizeData', 20);
plot_ellipse(Sigmax_T, xx_T, 'edgecolor', xx_color, 'confidence', 0.5);
plot_ellipse(Sigmax_T, xx_T, 'edgecolor', xx_color, 'confidence', 0.99);

% plot tangent line
tangent_color = hex2rgb('a64d79');

% plot vector from origin to xx_T
quiver(0, 0, xx_T(1), xx_T(2), 1, ...
    'color', tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);

x = linspace(-10, 10, 1000);
xx_T_N = xx_T./norm(xx_T);
y = 1/xx_T_N(2) - xx_T_N(1)/xx_T_N(2)*x;
plot(x, y, 'color',tangent_color); %light yellow

% plot points on Minkowski sum
scatter(xx_T_N(1), xx_T_N(2), 'MarkerFaceColor', tangent_color,...
 'MarkeredgeColor', tangent_color, 'SizeData', 30);

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
ax.XLim = [-2, 2];
ax.YLim = [-1.5, 2.5];
axis equal

%% lcc-center Zhu's equal: using bounding ellip without space transformation
figure; hold on;

% plot bounding ellip
plot_ellipse(Sigmaf, 'edgecolor', bounding_ellip_color);

% plot minksum
patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% Plot relative position error
% Plot relative position error
scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', xx_color, 'SizeData', 20);
plot_ellipse(Sigmax, xx, 'edgecolor', xx_color, 'confidence', 0.5);
plot_ellipse(Sigmax, xx, 'edgecolor', xx_color, 'confidence', 0.99);

% plot line from origin to xx
plot([0, xx(1)], [0, xx(2)], '--')

% plot intersection point
xx = s2.tc - s1.tc;
x_mink = xx./norm(Sigmaf^0.5 \ xx);
scatter(x_mink(1), x_mink(2), '*', 'SizeData', 150, 'MarkeredgeColor', hex2rgb('f44336'));

% plot tangent plane
tangent_color = hex2rgb('a64d79');
norm_plane = (Sigmaf \ xx) ./ norm(Sigmaf \ xx);
b_plane = x_mink'*norm_plane;

% plot norm at x_mink
quiver(x_mink(1), x_mink(2), norm_plane(1), norm_plane(2),5, ...
    'color', tangent_color, 'LineWidth', 0.8, 'MaxHeadSize', 1);

% plot the half space
x = -200:0.01:200;
y = b_plane/norm_plane(2) - norm_plane(1)/norm_plane(2)*x;
plot(x, y, 'color',tangent_color); %light yellow

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
ax.XLim = [-15, 15];
ax.YLim = [-12, 25];
axis equal

% %% lcc-center using closed-form Mink sum without space transformation
% plot line from origin to xx
plot([0, xx(1)], [0, xx(2)], '--')

% % Given objective, find x_mink(gradient, norm) 
% Constuct MinkSumClosedForm of s1 and s2
minkSum = MinkSumClosedForm(s1, s2, angle2rotm(s1.ang),angle2rotm(s2.ang));

% the intersection point should lie at the same line as origin to xx 
xx_norm = xx./norm(xx);

% % Optimization 1: use angular parameter as variable
% option = optimoptions('lsqnonlin',...
%             'Algorithm', 'levenberg-marquardt',...
%             'display', 'iter',...
%             'FunctionTolerance', 1e-8,...
%             'OptimalityTolerance', 1e-8);
% psi_opt = lsqnonlin(@(psi) func_lsq(psi, minkSum, xx_norm), 0,...
%     [], [], option);
% % Solution gradient in local frame
% m_opt = s1.GetGradientFromAngle(psi_opt);

% Optimization 2: use gradient as optimization variable
option = optimoptions('fmincon', 'Algorithm', 'interior-point',...
            'display', 'none');
m_opt = fmincon(@(m) func_con(m, minkSum, xx_norm), xx_norm,...
    [], [], [], [], [], [], @(m) nlcon(m, minkSum, xx_norm), option);

% Get x_mink
x_mink = minkSum.GetMinkSumFromGradient(m_opt);

scatter(x_mink(1), x_mink(2), '*', 'SizeData', 150, 'MarkeredgeColor', hex2rgb('f44336'));

% get m1 in current space
m1_current = angle2rotm(s1.ang)' \ m_opt;

% get norm
norm_plane = m1_current./norm(m1_current);

% plot tangent plane
new_tangent_color = hex2rgb('4faad1');

b_plane = x_mink'*norm_plane;

% plot norm at x_mink
quiver(x_mink(1), x_mink(2), norm_plane(1), norm_plane(2),5, ...
    'color', new_tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);

x = linspace(-40, 40, 1000);
y = b_plane/norm_plane(2) - norm_plane(1)/norm_plane(2)*x;
plot(x, y, 'color',new_tangent_color); %light yellow

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
ax.XLim = [-12, 12];
ax.YLim = [-12, 25];
axis equal

%% optimization functions
% Least squares optimization cost
function F = func_lsq(psi, minkSum, xx)
m = minkSum.s1.GetGradientFromAngle(psi);

x_mink = minkSum.GetMinkSumFromGradient(m);

F = norm(x_mink) - x_mink' * xx;
end

% Constrained optimization cost and constraint
% vector angle cost: find x_mink at the line xx
function F = func_con(m, minkObj, xx)
x_mink = minkObj.GetMinkSumFromGradient(m);

% Cost function
F = norm(x_mink) - x_mink' * xx;
end

% Nonlinear constraint: x_mink and xx angle less then pi/2
function [c,ceq] = nlcon(m, minkObj, xx)
x_mink = minkObj.GetMinkSumFromGradient(m);
% inequality constraint: x_mink and xx angle less then pi/2
% c = -1 * x_mink' * xx;
c=[];
% equality constraint: x on the surface
x = minkObj.s1.GetPointsFromGradient(m);
ceq = minkObj.s1.GetImplicitFunction(x);
end

