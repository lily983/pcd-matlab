function visualize_position_error(s1, s2, Sigma)
    n = 100;
    
    dimension = size(Sigma,1);
    
    xi_rand = mvnrnd(zeros(dimension,1), Sigma, n);
    
    if dimension == 2
        s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
    elseif dimension == 3
        s3 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});
    end

%     figure; hold on;
%     s3.PlotShape('r', 0.1);
%     s1.PlotShape('g', 0.3);

    for i = 1:size(xi_rand, 1)
        s3.tc = s2.tc + xi_rand(i,:)';
        s3.PlotShape(hex2rgb('ffd166')./255, 0.05);
    end
    
end