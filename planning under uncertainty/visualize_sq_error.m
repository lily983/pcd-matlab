% function visualize_sq_error(s1_param, s2_param, Sigma)
function visualize_sq_error(s2_param, Sigma)
    n = 10;
    
    dimension = size(Sigma,1);
    
    xi_rand = mvnrnd(zeros(dimension,1), Sigma, n);
    
%     s1 = SuperEllipse([s1_param.a(1), s1_param.a(2), s1_param.eps, 0, ...
%         s1_param.tc(1), s1_param.tc(2), s1_param.ang, 20]);
    
     s3 = SuperEllipse([s2_param.a(1), s2_param.a(2), s2_param.eps, 0, ...
        s2_param.tc(1), s2_param.tc(2), s2_param.ang, 20]);

%     figure; hold on;axis equal;axis off
    s3.PlotShape(hex2rgb('45AC59'), 0.4, 0.4);
%     s1.PlotShape('b', 0.3, 0.3);

    for i = 1:size(xi_rand, 1)
        s3.tc = s2_param.tc + xi_rand(i,:)';
        s3.PlotShape(hex2rgb('45AC59'), 0.1, 0.2);

    end
    
end