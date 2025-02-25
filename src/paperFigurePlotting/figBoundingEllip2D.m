close all;clear;clc;

% figure; hold on; axis equal; axis off;
figure; hold on; 

s1 = SuperEllipse([4, 10, 1, 0, 15, 18, 0.0, 50]);
s2 = SuperEllipse([5, 2, 1, 0, 18, 1, 0.3, 50]);

%%% Plot the effect of position errors
n = 20;

s1_color = hex2rgb('45AC59');
s1.PlotShape(s1_color, 1, 1);
Sigma1 = eye(2, 2) * 16 *1e-01;
x1_rand = mvnrnd(zeros(2,1), Sigma1, n);
s1_shift = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
    s1.tc(1), s1.tc(2), s1.ang, s1.N]);
for i = 1:size(x1_rand, 1)
    s1_shift.tc = s1.tc + x1_rand(i,:)';
    s1_shift.PlotShape(s1_color, 0.15, 0.0);
end

s2_color = hex2rgb('4FAAD1')
s2.PlotShape(s2_color, 1, 1);
Sigma2 =  eye(2, 2) * 8 *1e-01;
x2_rand = mvnrnd(zeros(2,1), Sigma2, n);
s2_shift = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
    s2.tc(1), s2.tc(2), s2.ang, s2.N]);
for i = 1:size(x1_rand, 1)
    s2_shift.tc = s2.tc + x2_rand(i,:)';
    s2_shift.PlotShape(s2_color, 0.15, 0.0);
end

%%% Visualize bounding ellipsoid
[mf, Sigmaf] = get_bounding_ellip(s1, s2);
[U, Lamda] = svd(Sigmaf);

s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
            s2.tc(1), s2.tc(2), s2.ang, s2.N]);
s3.a = sqrt(diag(Lamda));
s3.ang = rotm2angle(U);
s3.tc = mf;
s3.N=100;

% [~, Sigmaf_fp] = get_bounding_ellip_fixed_point(s1, s2);
% [U_fp, Lamda_fp] = svd(Sigmaf_fp);
% 
% s3_fp = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
%     s2.tc(1), s2.tc(2), s2.ang, s2.N]);
% s3_fp.a = sqrt(diag(Lamda_fp));
% s3_fp.ang = rotm2angle(U_fp);
% s3_fp.tc = mf;
% s3_fp.N=100;

% Get minkowski sum of s1 and s2
s1Points = s1.GetPoints()' - s1.tc';
s2Points = -1*s2.GetPoints()' + s2.tc';
pgon1 = polyshape(s1Points(:,1), s1Points(:,2));
pgon2= polyshape(s2Points(:,1), s2Points(:,2));
mSum = minkowskiSum(pgon1, pgon2);

patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('45498C'), 'FaceAlpha', 0.0, 'EdgeAlpha', 1);
s3.PlotShape(hex2rgb('ED5564'), 0.0, 0.8); % light red
% s3_fp.PlotShape(hex2rgb('ED5564'), 0.0, 1);

%%% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;

%% Shown Minkowski sum process
face_alpha = 0.4;
edge_alpha = 0.5;

s1_shift.tc = zeros(2,1);
s1_shift.PlotShape(s1_color, face_alpha, edge_alpha);

s2_shift.tc =[4;0];
s2_shift.PlotShape(s2_color, face_alpha, edge_alpha);

s2_shift.tc =[-4;0];
s2_shift.PlotShape(s2_color, face_alpha, edge_alpha);

s2_shift.tc = [-2.49396; -7.81831];
s2_shift.PlotShape(s2_color, face_alpha, edge_alpha);

s2_shift.tc = [2.49396; 7.81831];
s2_shift.PlotShape(s2_color, face_alpha, edge_alpha);
% 
% %%
% s1_gradient = s1.GetGradientsCanonical();
% minksum = MinkSumClosedForm(s1, s2, rot2(s1.ang), rot2(s2.ang));
% minksum_points = minksum.GetMinkSumFromGradient(s1_gradient);
% plot(minksum_points(1, :), minksum_points(2,:), 'color', hex2rgb('45498C'));

