% clc; clear; close all;

% Ellipsoids and position covariance matrix setting
s1 = SuperQuadrics({[0.6, 0.6, 1.2], [1,1], [0, 0]...
    [0;0;0], [1, 0,0,0], [20, 20]});
s2 = SuperQuadrics({[0.18, 0.18, 0.22],  [1,1], [0, 0]...
    [0.95; 0.95; 0], [1,0,0,0],[20,20]});

[flag, dist, ~, ~] = collision_cfc(s1,s2, 'constrained');

%     Covariance matrix of the position error distribution
Sigmax = zeros(3);
Sigmax(1,1) = 0.41;
Sigmax(2,2) = 0.41;
Sigmax(3,3) = 0.21;
% joint distribution
Sigmax = Sigmax.*2;

%     Get PCD values of each methods
PCDExact = getPCDR3(s1, s2, Sigmax, 'PCD-exact');
PCDConvex  =  getPCDR3(s1, s2, Sigmax, 'PCD-convex');
PCDEB = getPCDR3(s1, s2, Sigmax, 'PCD-EB-99');
PCDGMM  = getPCDR3(s1, s2, Sigmax, 'PCD-GMM');
PCDEllipBound = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-bound');
PCDEllipExact = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-exact');
PCDEllip2023Approximation = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-UB-approximation');
PCDEllip2023Exact = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-UB-exact');

% figure;hold on; axis equal
% s1.PlotShape('g', 0.5);
% s2.PlotShape('b', 0.5);
visualize_position_error(s1, s2, Sigmax);


