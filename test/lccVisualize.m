% This file visualize lcc related concepts in 2D
clc; clear; close all; 

%% lcc-center: ellipsoids case
figure; hold on; 

% Construct two ellipses
s1 = SuperEllipse([1, 6, 1, 0, -17, -12, 0.0, 50]);
s2 = SuperEllipse([5, 3, 1, 0, -12, 1, -0.5, 50]);

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

% plot the bounding ellipsoid for the Minkowksi sum
[mf, Sigmaf] = get_bounding_ellip(s1, s2);
[U_Sigmaf, Lamda_Sigmaf] = svd(Sigmaf);

s3 = SuperEllipse([s2.a(1), s2.a(2), 1, 0,...
            s2.tc(1), s2.tc(2), s2.ang, s2.N]);
s3.a = sqrt(diag(Lamda_Sigmaf));
s3.ang = rotm2angle(U_Sigmaf);
s3.tc = mf;
s3.N=100;
s3.PlotShape(hex2rgb('f44336'), 0.0, 1.0); %red

% plot position error
Sigma1 = [1, 0; 0, 9]*1e-01;
Sigma2 = [8, 0; 0, 8] * 1e-01;

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

% Here Sigmax_T is a diagonal matrix for simplification
sx_a = diag(sqrtm(Sigmax))';
sx = SuperEllipse([sx_a, 1, 0,...
            xx(1), xx(2), 0, s2.N]);

scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', hex2rgb('714bd2'), 'SizeData', 20);

% plot contour of position error
i=1;
while i<=4
    sx.PlotShape(hex2rgb('714bd2'), 0.0, 1.0); %purple
    i = i+1;
    sx.a = sx_a .* 1.5^i;
end

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
[U_T, Lamda_T] = svd(Sigmaf_T);
s3.a = sqrt(diag(Lamda_T));
s3.ang = rotm2angle(U_T);
% If transformed s3 is not a sphere, then check the transformation
s3.PlotShape(hex2rgb('f44336'), 0.0, 1.0); 

% transformed Minkowski sum
mSum_T=  Sigmaf^0.5 \ [mSum.Vertices(:,1), mSum.Vertices(:,2)]';
patch(mSum_T(1, :), mSum_T(2, :), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% transformed relative position error
xx_T = Sigmaf^0.5 \ xx;
scatter(xx_T(1), xx_T(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', hex2rgb('0b5394'), 'SizeData', 20);

Sigmax_T = Sigmaf^0.5 \ Sigmax / Sigmaf^0.5;

% Here Sigmax_T is not a diagonal matrix
[R_Sigmax_T, Lamda_Sigmax_T ,~] = svd(Sigmax_T);
sx_T_a = 1./sqrt(diag(Lamda_Sigmax_T))';

sx_T = SuperEllipse([sx_T_a, 1, 0,...
            xx_T(1), xx_T(2), rotm2angle(R_Sigmax_T), s2.N]);

% plot contour of position error
i=1;
while i<=3
    sx_T.PlotShape(hex2rgb('714bd2'), 0.0, 1.0); %purple
    i = i+1;
    sx_T.a =  sx_T_a .* 0.5^i;
end

% plot tangent line
tangent_color = hex2rgb('a64d79');

% plot vector from origin to xx_T
quiver(0, 0, xx_T(1), xx_T(2), 1, ...
    'color', tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);

x = linspace(-4, 6, 1000);
xx_T_N = xx_T./norm(xx_T);
y = 1/xx_T_N(2) - xx_T_N(1)/xx_T_N(2)*x;
plot(x, y, 'color',tangent_color); %light yellow

% plot points on Minkowski sum
scatter(xx_T_N(1), xx_T_N(2), 'MarkerFaceColor', tangent_color,...
 'MarkeredgeColor', tangent_color, 'SizeData', 30);

% plot the half space
[xymi,xymax]=bounds([x; y]',1);
x_h=linspace(xymi(1),xymax(1),33);
y_h=linspace(xymi(2),xymax(2),33);
[X_H,Y_H]=meshgrid(x_h, y_h);
XY=[X_H(:),Y_H(:)].';
inside = all(xx_T_N' * XY -1<=0, 1);
plot(XY(1,inside),XY(2,inside),'.', 'color', tangent_color);

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
ax.XLim = [-4, 6];
ax.YLim = [-3, 6];
axis equal

%% lcc-center using closed-form Mink sum without space transformation
figure; hold on;

% plot bounding ellip
s3.a = sqrt(diag(Lamda_Sigmaf));
s3.ang = rotm2angle(U_Sigmaf);
s3.PlotShape(hex2rgb('f44336'), 0.0, 1.0); 

% plot minksum
patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% Plot relative position error
scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', hex2rgb('714bd2'), 'SizeData', 20);

sx = SuperEllipse([1./diag(sqrtm(Sigmax))', 1, 0,...
            xx(1), xx(2), 0, s2.N]);
% plot contour of position error
i=1;
while i<=4
    sx.PlotShape(hex2rgb('714bd2'), 0.0, 1.0); %purple
    i = i+1;
    sx.a = sx_a .* 1.5^i;
end

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
quiver(x_mink(1), x_mink(2), norm_plane(1) + x_mink(1), norm_plane(2) + x_mink(2),1, ...
    'color', tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);

x = linspace(-8, 20, 1000);
y = b_plane/norm_plane(2) - norm_plane(1)/norm_plane(2)*x;
plot(x, y, 'color',tangent_color); %light yellow

% plot the half space
[xymi,xymax]=bounds([x; y]',1);
x_h=linspace(xymi(1),xymax(1),40);
y_h=linspace(xymi(2),xymax(2),40);
[X_H,Y_H]=meshgrid(x_h, y_h);
XY=[X_H(:),Y_H(:)].';
inside = all(norm_plane' * XY - b_plane<=0, 1);
plot(XY(1,inside),XY(2,inside),'.', 'color', tangent_color);

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
ax.XLim = [-8, 11];
ax.YLim = [-10, 18];
axis equal