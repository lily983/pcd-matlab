function visualize_position_error(s1, s2, Sigma)
    n = 40;
    
    dimension = size(Sigma,1);
    
    xi_rand = mvnrnd(zeros(dimension,1), Sigma, n);
    
    if dimension == 2
        s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
    elseif dimension == 3
        s3 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});
    end

    figure; hold on;axis equal;axis off
    s3.PlotShape('r', 0.4, 0.4);
    s1.PlotShape('b', 0.3, 0.3);

    for i = 1:size(xi_rand, 1)
        s3.tc = s2.tc + xi_rand(i,:)';
        s3.PlotShape(hex2rgb('45AC59'), 0.2, 0.4);
%         pause(0.5)
    end
    
end