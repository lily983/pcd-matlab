function visualize_3D_case(result1, result2, col, Cov1, Cov2)

N=[20, 20];

s1 = SuperQuadrics({result1(1:3,col)', result1(11:12,col)', [0, 0]...
result1(8:10,col), result1(4:7,col)', N});

s2 = SuperQuadrics({result2(1:3,col)', result2(11:12,col)', [0, 0]...
   result2(8:10,col), result2(4:7,col)', N});

% visualize_bounding_ellip(s1, s2);

visualizePositionErrors(s1, s2, Cov1, Cov2)

% Transform position error covariance matrix to world space
R1 = quat2rotm(s1.q);
Cov1 = R1 * Cov1 * R1';
R2 = quat2rotm(s2.q);
Cov2 = R2 * Cov2 * R2';

% relative position error
xx = s2.tc-s1.tc;
Covx = Cov1 + Cov2;

% [p_center, ~, a_center, x_center] = linearChanceConstraintSQ(s1, s2, xx, Covx, 'center-point', false)

% [p_center_cfc, ~, a_center_cfc, x_center_cfc] = linearChanceConstraintSQ(s1, s2, xx, Covx, 'center-point-cfc', false)

% [p_tangent_cfc, ~, a_tangent_cfc, x_tangent_cfc] = linearChanceConstraintSQ(s1, s2, xx, Covx, 'tangent-point-cfc', false)

% [prob_exact, ~] = exactProbTranslation(s1, s2,  Cov2, 1e+03);
 
[prob_diverg, ~] = divergenceMesh(s1, s2, xx, Covx)

% [prob, time] = exactProbTranslationTwoErrors(s1, s2, Cov1, Cov2, 1e+03)

%  [prob, time] = fastExactProbTranslationTwoErrors(s1, s2, Cov1, Cov2, 1e+04)

% %% This figure plots s1, s2 and their bounding ellip
% figure; hold on
% s1.PlotShape('b', 0.7);
% s2.PlotShape('g', 0.7);
% 
% % Get bounding ellipsoid for each object
% s1Points = s1.GetPoints();
% s2Points = s2.GetPoints();
% 
% [invCov1, c1] = MinVolEllipse(s1Points, 0.001);
% [invCov2, c2] = MinVolEllipse(s2Points, 0.001);
% Cov1 = inv(invCov1);
% Cov2 = inv(invCov2);
% 
% plot_ellipse(Cov1, c1, 'edgecolor', 'r');
% plot_ellipse(Cov2, c2, 'edgecolor', 'r');
% 
% This figure plots the exact minkowski sum, bounding ellip, tangent
% plane, position error ellip
figure; hold on

% Get bounding ellip for the Minkowski sum of two ellips
% Sigmaf = (1 + sqrt(trace(Cov2)/trace(Cov1))) * Cov1...
%     + (1 + sqrt(trace(Cov1)/trace(Cov2))) * Cov2;
% 
% % Plot bounding ellip
% plot_ellipse(Sigmaf, xx, 'confidence', 0.99, 'edgecolor', 'r');

% Plot Minkowski Sum
m1 = s1.GetGradientsCanonical();
minkSum = MinkSumClosedForm(s1, s2, quat2rotm(s1.q), quat2rotm(s2.q));
x_mink = minkSum.GetMinkSumFromGradient(m1);
N = s1.N;
X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

mink_sf = surf(X, Y, Z,...
 'EdgeColor', '#93c47d', 'EdgeAlpha', 0.0,...
 'FaceAlpha', 0.8, 'FaceColor', 'b'); %plot contact space

% Plot tangent plane
% plotPlane(a_center, x_center, 'b')
% plotPlane(a_center_cfc, x_center_cfc, 'm')
plotPlane(a_tangent_cfc, x_tangent_cfc, 'g')

% 
% 
% plot position error ellipse
position_error_color = hex2rgb('714bd2');
scatter3(xx(1), xx(2), xx(3), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);
% plot_ellipse(inv(Covx), xx, 'confidence', 0.99, 'edgecolor', position_error_color)
% 
% %% This figure plots minkowski sum, bounding ellip, and position error
% % ellip.
% figure; hold on; axis equal
% % Plot bounding ellip
% plot_ellipse(Sigmaf, 'confidence', 0.99, 'edgecolor', 'r');
% 
% mink_sf = surf(X, Y, Z,...
%  'EdgeColor', '#93c47d', 'EdgeAlpha', 0.0,...
%  'FaceAlpha', 0.8, 'FaceColor', 'b'); %plot contact space
% 
% % plot position error ellipse
% plot_ellipse(Covx, xx, 'confidence', 0.99, 'edgecolor', position_error_color)

end