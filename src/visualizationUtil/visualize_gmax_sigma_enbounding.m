function visualize_gmax_sigma_enbounding(s1, s2, result, col, Sigma_g)
    
    mu_vee = result(29:34, col);
    g_max_vee = result(22:27, col);
    mu = get_SE3_matrix(mu_vee);
    g_max = get_SE3_matrix(g_max_vee);
    
    s2.tc = mu(1:3,4);
    s2.q = rotm2quat(mu(1:3,1:3));

    s3 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});
    
    close all;
    % Visualize sigma effect
    figure; hold on;
    s1.PlotShape('b', 0.3)
    s2.PlotShape('g', 0.3)
    pause(1)
    
    n = 100;
    gi_vee_rand = mvnrnd(get_vee_vector(mu), Sigma_g, n);

    for i=1:size(gi_vee_rand, 1)
        gi_matrix = get_SE3_matrix(gi_vee_rand(i,1:6)');
        gi_tc = gi_matrix(1:3,4);
        gi_q = rotm2quat(gi_matrix(1:3,1:3));
        s3.tc = gi_tc;
        s3.q = gi_q;
        s3.PlotShape('g', 0.1);
    end
    pause(1)
    
    % Visualize gmax
    figure; hold on;
    s1.PlotShape('b', 0.3)
    s2.PlotShape('g', 0.3)
    
    s3.tc = g_max(1:3,4);
    s3.q = rotm2quat(g_max(1:3,1:3));
    s3.PlotShape('r',0.1);
    pause(1)
    
    % Visualize bounding superquadrics
    figure; hold on;
    s1.PlotShape('b', 0.3)
    s2.PlotShape('g', 0.3)
    sq_patch = enlarged_bounding_superquadrics(s2, mu, Sigma_g);
    
    
    


end