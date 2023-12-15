function visualizePositionErrors(s1, s2, Sigma1, Sigma2)
    n = 20;
    
    dimension = size(Sigma1,1);
    
    x1_rand = mvnrnd(zeros(dimension,1), Sigma1, n);
    
    x2_rand = mvnrnd(zeros(dimension,1), Sigma2, n);
    
    if dimension == 2
        s3 = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
        s1.tc(1), s1.tc(2), s1.ang, s1.N]);
        s4 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
    elseif dimension == 3
        s3 = SuperQuadrics({s1.a, s1.eps, s1.taper, s1.tc, s1.q, s1.N});
        s4 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});
    end

    figure; hold on;axis equal;axis off
    s3.PlotShape('r', 0.4);
    s1.PlotShape('b', 0.3);

    for i = 1:size(x1_rand, 1)
        s3.tc = s1.tc + x1_rand(i,:)';
        s3.PlotShape(hex2rgb('9fc5e8'), 0.1);
        
        s4.tc = s2.tc + x2_rand(i,:)';
        s4.PlotShape(hex2rgb('45AC59'), 0.1);
%         pause(0.5)
    end
    
end