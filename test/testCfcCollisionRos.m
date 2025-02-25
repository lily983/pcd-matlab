close all; clc; clear;

a1=[0.06723, 0.080807, 0.23949];
pos1=[0.468383, -0.245903, 0.42022];
q1=[-0.0767502, 0.847322, 0.0153524, 0.52528];

a2=[0.12735, 0.16917, 0.19533];
pos2=[0.69906, -1.37e-05, 0.16];
q2=[0.999998, 5.57677e-05, -0.000246612, 0.0018466];


SN = [10, 10];
s1 = SuperQuadrics({a1, [1,1], [0, 0]...
    pos1', q1, SN});
s2 = SuperQuadrics({a2, [1,1], [0, 0]...
    pos2', q2, SN});

cov1=zeros(3);
cov2=eye(3);
cov2(1,1) = 0.00021;
cov2(2,2) = 0.00024;
cov2(3,3) = 0.00027;
% visualizePositionErrors(s1, s2, cov1, cov2)

mx = (pos2-pos1)';
Sigmax = cov1+cov2;

[prob_center, ~] = linearChanceConstraintSQ(s1, s2, mx, Sigmax, 'center-point', false)

[prob_tangent, ~] = linearChanceConstraintSQ(s1, s2, mx, Sigmax, 'tangent-point-cfc', false)

[prob_exact, ~] = exactProbTranslation(s1, s2, cov2, 1e+03)

[~, dist, pt_cls, ~] = collision_cfc(s1, s2, 'least-squares')

%%
figure; hold on;axis equal;axis off
% s1.PlotShape(hex2rgb('9fc5e8'), 0.6);
% s2.PlotShape(hex2rgb('45AC59'), 0.8);

m1 = s1.GetGradientsCanonical();
minkSum = MinkSumClosedForm(s1, s2,quat2rotm(s1.q),quat2rotm(s2.q));
exact_mink = minkSum.GetMinkSumFromGradient(m1);

% plot the exact minksum
N = s1.N;
X = reshape(exact_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(exact_mink(2,:), N(1), N(2));
Z = reshape(exact_mink(3,:), N(1), N(2));

mink_sf = surf(X, Y, Z,...
 'EdgeColor', '#93c47d', 'EdgeAlpha', 1,...
 'FaceAlpha', 0.8, 'FaceColor', '#93c47d'); %plot contact space

position_error_color = hex2rgb('714bd2');
scatter3(mx(1), mx(2), mx(3), 'MarkerFaceColor', position_error_color,...
 'MarkeredgeColor', position_error_color, 'SizeData', 20);

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% ax.ZAxisLocation = 'origin';
ax.FontSize = 8;
% axis equal
