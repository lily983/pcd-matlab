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
[invSigmaf, mf] = MinVolEllipse(mSum.Vertices', 0.001)
Sigmaf = inv(invSigmaf);
[R_Sigmaf, Lamda_Sigmaf,~] = svd(Sigmaf);

s3_a = sqrt(diag(Lamda_Sigmaf))';
s3 = SuperEllipse([s3_a, 1, 0,...
            mf(1), mf(2), rotm2angle(R_Sigmaf), s2.N]);
        
% plot bounding ellipse
s3.PlotShape(hex2rgb('f44336'), 0.0, 1.0); 

% plot position error
Sigma1 = [4, 0; 0, 9]*1e-01;
Sigma2 = [8, 0; 0, 15] * 1e-01;

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

% Here Sigmax is simplified to be a diagonal matrix, so we can use diag()
% to get the value
[R_Sigmax, Lamda_Sigmax,~] = svd(Sigmax);
sx_a = sqrt(diag(Lamda_Sigmax))';

sx = SuperEllipse([sx_a, 1, 0,...
            xx(1), xx(2), rotm2angle(R_Sigmax), s2.N]);

scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', hex2rgb('714bd2'), 'SizeData', 20);
i=1;
while i<=3
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

% First get an bounding ellipsoid for the Minkowski sum
[invSigmaf, mf] = MinVolEllipse(mSum.Vertices', 0.001)
Sigmaf = inv(invSigmaf);
[R_Sigmaf, Lamda_Sigmaf,~] = svd(Sigmaf);

s3_a = sqrt(diag(Lamda_Sigmaf))';
s3 = SuperEllipse([s3_a, 1, 0,...
            mf(1), mf(2), rotm2angle(R_Sigmaf), s2.N]);
        
% transformed bounding ellip
Sigmaf_T = Sigmaf^0.5 \ Sigmaf / Sigmaf^0.5;
[R_Sigmaf_T, Lamda_Sigmaf_T] = svd(Sigmaf_T);

s3.a = sqrt(diag(Lamda_Sigmaf_T));
s3.ang = rotm2angle(R_Sigmaf_T);
s3.tc = Sigmaf^0.5 \ s3.tc;
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

% Here Sigmax is not a diagonal matrix
[R_Sigmax_T, Lamda_Sigmax_T ,~] = svd(Sigmax_T);
sx_T_a = 1./sqrt(diag(Lamda_Sigmax_T))';
sx_T = SuperEllipse([sx_T_a, 1, 0,...
            xx_T(1), xx_T(2), rotm2angle(R_Sigmax_T), s2.N]);
        
% plot contour relative position error
i=6;
while i<=8
    i = i+1;
    sx_T.a =  sx_T_a .* 0.7^i;
    sx_T.PlotShape(hex2rgb('0b5394'), 0.0, 1.0); %purple
end

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

% plot the half space
[xymi,xymax]=bounds([x; y]',1);
x_h=linspace(xymi(1),xymax(1),50);
y_h=linspace(xymi(2),xymax(2),50);
[X_H,Y_H]=meshgrid(x_h, y_h);
XY=[X_H(:),Y_H(:)].';
inside = all(xx_T_N' * XY -1<=0, 1);
plot(XY(1,inside),XY(2,inside),'.', 'color', tangent_color);

% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;
ax.XLim = [-1.5, 1.5];
ax.YLim = [-1.5, 2.5];
axis equal

%% lcc-center Zhu's equal: using bounding ellip without space transformation
figure; hold on;

% plot bounding ellip
s3.a = s3_a;
s3.ang = rotm2angle(R_Sigmaf);
s3.PlotShape(hex2rgb('f44336'), 0.0, 1.0); 

% plot minksum
patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% Plot relative position error
scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', hex2rgb('714bd2'), 'SizeData', 20);

i=1;
while i<=3
    sx.PlotShape(hex2rgb('714bd2'), 0.0, 1.0); %purple
    i = i+1;
    sx.a = sx_a .* 1.5^i;
end

% plot line from origin to xx
plot([0, xx(1)], [0, xx(2)], '--')

% plot intersection point
% here x_mink is by observation
x_mink = [3.2; 8.97357];
scatter(x_mink(1), x_mink(2), '*', 'SizeData', 150, 'MarkeredgeColor', hex2rgb('f44336'));

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

x = linspace(-200, 200, 5000);
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
ax.XLim = [-12, 12];
ax.YLim = [-12, 25];
axis equal

%% lcc-center using closed-form Mink sum without space transformation
figure; hold on;

% plot bounding ellip
s3.a = s3_a;
s3.ang = rotm2angle(R_Sigmaf);
s3.PlotShape(hex2rgb('f44336'), 0.0, 1.0); 

% plot minksum
patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('93c47d'), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.0);

% Plot relative position error
scatter(xx(1), xx(2), 'MarkerFaceColor', hex2rgb('714bd2'),...
 'MarkeredgeColor', hex2rgb('714bd2'), 'SizeData', 20);

i=1;
while i<=3
    sx.PlotShape(hex2rgb('714bd2'), 0.0, 1.0); %purple
    i = i+1;
    sx.a = sx_a .* 1.5^i;
end

% plot line from origin to xx
plot([0, xx(1)], [0, xx(2)], '--')

% plot intersection point
% here x_mink is by observation
x_mink = [3.2; 8.97357];
scatter(x_mink(1), x_mink(2), '*', 'SizeData', 150, 'MarkeredgeColor', hex2rgb('f44336'));

% get norm at x_mink
norm_plane

% plot tangent plane
tangent_color = hex2rgb('a64d79');

b_plane = x_mink'*norm_plane;

% plot norm at x_mink
quiver(x_mink(1), x_mink(2), norm_plane(1) + x_mink(1), norm_plane(2) + x_mink(2),1, ...
    'color', tangent_color, 'LineWidth', 1, 'MaxHeadSize', 1);

x = linspace(-40, 40, 1000);
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
ax.XLim = [-12, 12];
ax.YLim = [-12, 25];
axis equal