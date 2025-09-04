function visualize_2D_case(result, col, Sigma)
% Visualize col-th case in result
% Use this function with test_sphere_2D.m
% Type: sphere or ellipsoid
s1 = SuperEllipse([result(6, col), result(6, col), 1, 0 ...
        result(7, col), result(8, col), result(9, col), 25]);
    
s2 = SuperEllipse([result(11, col), result(11, col), 1, 0 ...
        result(12, col), result(13, col), result(14, col), 25]);
    
    % visualize bounding ellip and minkowski sum
    visualize_bounding_ellip(s1, s2);
    visualize_position_error(s1, s2, Sigma);
%     flagNow = collision_mesh(s1, s2);
    
    % visualize gaussian mixture function value
    [mf, Sigmaf] = get_bounding_ellip(s1, s2);
%     visualize_gmm_integration_2D(mf, Sigmaf, Sigma)
    visualize_gmm_2D(mf, Sigmaf, 10,10,3)
    
    % visualize error pdf distribution
%     visualize_error_pdf([0;0], Sigma);
    
end