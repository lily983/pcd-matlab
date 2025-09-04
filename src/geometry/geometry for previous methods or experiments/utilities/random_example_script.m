%% Generate superquadric and sampling points
clear
close all

% generate a random superquadric point cloud
% generating superquadric parameters
epsilon = 2 * rand(1, 2) + 0.004;
a = 2 * rand(1, 3) + 0.5;
euler = [2 * pi * rand, 2 * pi * rand, pi * rand] ;
t = 0.2 * rand(1, 3) - 0.1;
R = eul2rotm(euler);
x_gt = [epsilon, a, euler, t];

% point cloud sampling arclength
arclength = 0.2;

% generating random partial view
partial_ratio = 0.4;
point = randomPartialSuperquadrics(x_gt, arclength, partial_ratio);

% add noise
noise_level = 0.1;
noise = rand(3, size(point, 2)) * a(1) * noise_level - a(1) * noise_level / 2;
point = point + noise;

% add outlier
outlier_ratio = 0.1;
outlier = mvnrnd(t, 2 .* eye(3), floor(outlier_ratio * size(point, 2)))';
point = [point, outlier];

figure(1)
showPoints(point)
axis equal
title('Partial Point Cloud with Noise and Outliers')

%% Superquadric Recovery

[x_ems] = EMS(point, 'OutlierRatio', 0.2, 'IncludeAllInliers', true, 'InlierThreshold', 0.9);

disp('---------------------------------------------------------------------')
disp('Groud Truth parameters are');
disp(x_gt)
disp('---------------------------------------------------------------------')
disp('EMS Fitted parameters are')
disp(x_ems)
disp('---------------------------------------------------------------------')

figure(2)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_ems, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.1, 'Light', 1);
hold off
title('EMS Recovered Superquadric')