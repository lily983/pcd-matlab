function visualizePositionErrors(s1, s2, Cov1, Cov2, color1, color2)
if nargin == 4
    color1 = hex2rgb('45AC59');
    color2= hex2rgb('9fc5e8');
end
    n = 7;
    
    dimension = size(Cov1,1);
    
    % Transform position error covariance matrix to world space
    R1 = quat2rotm(s1.q);
    Cov1 = R1 * Cov1 * R1';
    R2 = quat2rotm(s2.q);
    Cov2 = R2 * Cov2 * R2';
    
    x1_rand = mvnrnd(zeros(dimension,1), Cov1, n);
    
    x2_rand = mvnrnd(zeros(dimension,1), Cov2, n);
    
    if dimension == 2
        s3 = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
        s1.tc(1), s1.tc(2), s1.ang, s1.N]);
        s4 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
    elseif dimension == 3
        s3 = SuperQuadrics({s1.a, s1.eps, s1.taper, s1.tc, s1.q, s1.N});
        s4 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});
    end

%     figure; hold on;axis equal;axis off
%     s1.PlotShape(color1, 0.6);
%     s2.PlotShape(color2, 0.8);

    for i = 1:size(x1_rand, 1)
        s3.tc = s1.tc + x1_rand(i,:)';
        s3.PlotShape(color1, 0.1, 0.0);
        
        s4.tc = s2.tc + x2_rand(i,:)';
        s4.PlotShape(color2, 0.1, 0.0);
%         pause(0.5)
    end
    
end