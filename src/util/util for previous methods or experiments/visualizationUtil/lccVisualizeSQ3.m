% This file visualize lcc related concepts in 2D
% This file uses bounding ellipsoid to bound s1 and s2, then use liu2023 to get lcc-center
clc; clear; close all; 

%% lcc-center: SQ case
figure; hold on; 

% Construct two SQ
s1 = SuperEllipse([2, 6, 0.2, 0, -19, -9, -0.5, 100]);
s2 = SuperEllipse([5, 3, 0.2, 0, -14, 5, -0.0, 100]);

s1_color = hex2rgb('45AC59'); %green
s1.PlotShape(s1_color, 1);

s2_color = hex2rgb('4FAAD1'); %blue
s2.PlotShape(s2_color, 1);

% Get bounding ellipsoid for s1 and s2
s1Points = s1.GetPoints();
s2Points = s2.GetPoints();
[invSigma1, e1] = MinVolEllipse(s1Points, 0.001);
[invSigma2, e2] = MinVolEllipse(s2Points, 0.001);
Sigma1 = inv(invSigma1);
Sigma2 = inv(invSigma2);

% plot bounding ellips for s1 and s2
plot_ellipse(Sigma1, e1, 'edgecolor', s1_color);
plot_ellipse(Sigma2, e2, 'edgecolor', s2_color);

% Get bounding ellipse for the Minkowski sum of the two ellipses
Sigmaf = (1 + sqrt(trace(Sigma2)/trace(Sigma1))) * Sigma1 + (1 + sqrt(trace(Sigma1)/trace(Sigma2))) * Sigma2;

% Plot bounding ellipse
bounding_ellip_color = hex2rgb('f44336');
plot_ellipse(Sigmaf, 'edgecolor', bounding_ellip_color);

% Get the exact Minkowski sum of two SQs
m1 = s1.GetGradientsCanonical();
minkSum = MinkSumClosedForm(s1, s2, angle2rotm(s1.ang), angle2rotm(s2.ang));
exact_mink = minkSum.GetMinkSumFromGradient(m1);

% plot the exact minksum
exact_mink_color = hex2rgb('93c47d');
patch(exact_mink(1,:), exact_mink(2,:), exact_mink_color, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% plot position error
% here we define the position error covariance is defined in the world
% frame, so we apply rotation
Sigmax1 = angle2rotm(s1.ang) * [4, 0; 0, 9]*1e-01 * angle2rotm(s1.ang)';
Sigmax2 = angle2rotm(s2.ang) * [8, 0; 0, 15] * 1e-01 * angle2rotm(s2.ang)';

n = 20;
x1_rand = mvnrnd(zeros(2,1), Sigmax1, n);
s1_shift = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
    s1.tc(1), s1.tc(2), s1.ang, s1.N]);
for i = 1:size(x1_rand, 1)
    s1_shift.tc = s1.tc + x1_rand(i,:)';
    s1_shift.PlotShape(s1_color, 0.15, 0.0);
end

x2_rand = mvnrnd(zeros(2,1), Sigmax2, n);
s2_shift = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
    s2.tc(1), s2.tc(2), s2.ang, s2.N]);
for i = 1:size(x1_rand, 1)
    s2_shift.tc = s2.tc + x2_rand(i,:)';
    s2_shift.PlotShape(s2_color, 0.15, 0.0);
end

% Plot the relative position error
xx = s2.tc - s1.tc;
Sigmax =Sigmax1+Sigmax2;

% plot position error ellipse
position_error_color = hex2rgb('714bd2');
scatter(xx(1), xx(2), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);
plot_ellipse(Sigmax, xx, 'confidence', 0.5, 'edgecolor', position_error_color)
plot_ellipse(Sigmax, xx, 'confidence', 0.99, 'edgecolor', position_error_color)

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
axis equal

%% lcc-center: Zhu's method
% % space transformation inv(Sigmaf)^0.5
% figure; hold on;
%         
% % transformed bounding ellip
% Sigmaf_T = Sigmaf^0.5 \ Sigmaf / Sigmaf^0.5;
% [R_Sigmaf_T, Lamda_Sigmaf_T] = svd(Sigmaf_T);
% 
% % plot transformed bounding ellip(should be a sphere)
% plot_ellipse(Sigmaf_T, 'edgecolor', bounding_ellip_color);
% 
% % transformed exact Minkowski sum
% mSum_T=  Sigmaf^0.5 \ [exact_mink(1,:); exact_mink(2,:)];
% patch(mSum_T(1, :), mSum_T(2, :), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);
% 
% % transformed relative position error
% xx_T = Sigmaf^0.5 \ xx;
% Sigmax_T = Sigmaf^0.5 \ Sigmax / Sigmaf^0.5;
% 
% % plot transformed position error
% scatter(xx_T(1), xx_T(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
%  'MarkeredgeColor', position_error_color, 'SizeData', 20);
% plot_ellipse(Sigmax_T, xx_T, 'edgecolor', position_error_color, 'confidence', 0.5);
% plot_ellipse(Sigmax_T, xx_T, 'edgecolor', position_error_color, 'confidence', 0.99);
% 
% % plot tangent line
% tangent_color = hex2rgb('a64d79');
% 
% % plot vector from origin to xx_T
% quiver(0, 0, xx_T(1), xx_T(2), 1, ...
%     'color', tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);
% 
% x = linspace(-10, 10, 1000);
% xx_T_N = xx_T./norm(xx_T);
% y = 1/xx_T_N(2) - xx_T_N(1)/xx_T_N(2)*x;
% plot(x, y, 'color',tangent_color); %light yellow
% 
% % plot points on Minkowski sum
% scatter(xx_T_N(1), xx_T_N(2), 'MarkerFaceColor', tangent_color,...
%  'MarkeredgeColor', tangent_color, 'SizeData', 30);
% 
% % Put axes center at the origin
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% ax.FontSize = 8;
% ax.XLim = [-2, 2];
% ax.YLim = [-1.5, 2.5];
% axis equal

%% lcc-center Zhu's equal: using bounding ellip without space transformation
figure; hold on;

% plot bounding ellip
plot_ellipse(Sigmaf, 'edgecolor', bounding_ellip_color);

% plot minksum
patch(exact_mink(1,:), exact_mink(2,:), exact_mink_color, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% Plot relative position error
scatter(xx(1), xx(2), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);
plot_ellipse(Sigmax, xx, 'edgecolor', position_error_color, 'confidence', 0.5);
plot_ellipse(Sigmax, xx, 'edgecolor', position_error_color, 'confidence', 0.99);

% plot line from origin to xx
plot([0, xx(1)], [0, xx(2)], '--')

% plot intersection point
x_mink = xx./norm(Sigmaf^0.5 \ xx);
scatter(x_mink(1), x_mink(2), '*', 'SizeData', 150, 'MarkeredgeColor', position_error_color);

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

% % Given objective, find x_mink(gradient, norm)
% the intersection point should lie at the same line as origin to xx 
xx_norm = xx./norm(xx);

% Optimization 1: optimize on angular parameters

% option = optimoptions('lsqnonlin',...
%             'Algorithm', 'levenberg-marquardt',...
%             'display', 'iter',...
%             'FunctionTolerance', 1e-8,...
%             'OptimalityTolerance', 1e-8);
% 
% psi_opt = lsqnonlin(@(psi) func_lsq(psi, minkSum, xx_norm), 0,...
%     [], [], option);
% % Solution gradient in local frame
% m_opt = s1.GetGradientFromAngle(psi_opt);

% Optimization 2: optimize on gradient. constraint is that x(m) is on
% the Minkowski sum boundary

option = optimoptions('fmincon', 'Algorithm', 'interior-point',...
            'display', 'none');
m_opt = fmincon(@(m) func_con(m, minkSum, xx_norm), xx_norm,...
    [], [], [], [], [], [], @(m) nlcon(m, minkSum, xx_norm), option);

% Get x_mink
x_mink = minkSum.GetMinkSumFromGradient(m_opt);

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

scatter(x_mink(1), x_mink(2), '*', 'SizeData', 150, 'MarkeredgeColor', new_tangent_color);

% lcc-tangent-cfm
xx_T = Sigmax^0.5 \ xx;

% Get the transformed exact Minkowski sum of two SQs
m1 = s1.GetGradientsCanonical();
minkSum_T = MinkSumClosedForm(s1, s2, Sigmax^0.5 \ angle2rotm(s1.ang), Sigmax^0.5 \ angle2rotm(s2.ang));
exact_mink_T = minkSum_T.GetMinkSumFromGradient(m1);

option = optimoptions('lsqnonlin',...
                    'Algorithm', 'levenberg-marquardt',...
                    'display', 'none',...
                    'FunctionTolerance', 1e-8,...
                    'OptimalityTolerance', 1e-8);
psi_opt = lsqnonlin(@(psi) func_lsq_tangent(psi, minkSum_T, xx_T), 0,...
            [], [], option);
        
m_opt = s1.GetGradientFromAngle(psi_opt);

x_mink_T = minkSum_T.GetMinkSumFromGradient(m_opt);

m_tangent = angle2rotm(s1.ang)' \ m_opt; 

a_tangent = m_tangent ./ norm(m_tangent);

x_mink_tangent = Sigmax^0.5 * x_mink_T;

b_tangent = x_mink_tangent' * a_tangent;

% lcc-tangent-cfm color
lcc_tangent_color = hex2rgb('1d5406');

% plot norm at x_mink
quiver(x_mink_tangent(1), x_mink_tangent(2), a_tangent(1), a_tangent(2), 5, ...
    'color', lcc_tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);

x = linspace(-40, 40, 1000);
y = b_tangent/a_tangent(2) - a_tangent(1)/a_tangent(2)*x;
plot(x, y, 'color',lcc_tangent_color); %light yellow

scatter(x_mink_tangent(1), x_mink_tangent(2), '*', 'SizeData', 150, 'MarkeredgeColor', lcc_tangent_color);

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

% Least squares optimization cost for lcc-tangent
function F = func_lsq_tangent(psi, minkSum_T, xx_T)
m = minkSum_T.s1.GetGradientFromAngle(psi);

xMink = minkSum_T.GetMinkSumFromGradient(m);

F = 0.5 * sum((xx_T - xMink).^2);
end