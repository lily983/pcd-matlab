% Let s1 be robot, s2 be obstacle
% The distance between s1 and s2 and the covariance matrix varies
s1 = SuperQuadrics({[0.3,0.3,0.3], [1,1], [0, 0]...
    [0;0;0], [1, 0,0,0], [20, 20]});

center_A = 0.8.*ones(3,1);
covariance_A = 0.1.*eye(3);
s2 = SuperQuadrics({[0.5,0.5,0.5],  [1,1], [0, 0]...
    center_A, [1,0,0,0],[20,20]});
PCDExact = getPCDR3(s1, s2, covariance_A, 'PCD-exact');
PCDConvex  =  getPCDR3(s1, s2, covariance_A, 'PCD-convex');
PCDEB = getPCDR3(s1, s2, covariance_A, 'PCD-EB-99');
PCDGMM  = getPCDR3(s1, s2, covariance_A, 'PCD-GMM');
PCDEllipBound = getPCDR3(s1, s2, covariance_A, 'PCD-ellip-bound');
PCDEllipExact = getPCDR3(s1, s2, covariance_A, 'PCD-ellip-exact');
PCDEllip2023Approximation = getPCDR3(s1, s2, covariance_A, 'PCD-ellip-UB-approximation');
PCDEllip2023Exact = getPCDR3(s1, s2, covariance_A, 'PCD-ellip-UB-exact');
visualize_position_error(s1, s2, covariance_A);
