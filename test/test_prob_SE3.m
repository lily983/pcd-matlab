%% Test 1: SE(3) pose error case. prob_g_max, prob_g_s1, exact_prob_g, exact_prob_x, prob_x_poly, prob_x_s1
close all;clear;clc;
add_path()

N = [20,20];

s1 = SuperQuadrics({0.02*ones(1,3), [0.8,0.7], [0, 0]...
zeros(3,1), [1, 0, 0, 0], N});

for i=1:10
   
    init_quat = rand(1,4);
    init_quat = init_quat/norm(init_quat);
    init_posi = 0.09*rand(3,1);
    s2 = SuperQuadrics({[0.01,0.02,0.03], [0.4,0.8], [0, 0]...
        init_posi, init_quat, N});

    [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'fixed-point');
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =10.0000e-06;
    Sigma_g(2,2) =9.0000e-06;
    Sigma_g(3,3) =8.0000e-06;
    Sigma_g(4,4) =1.0000e-06;
    Sigma_g(5,5) =4.0000e-06;
    Sigma_g(6,6) =9.0000e-06;
   
    % Compute prob_g_max
    [pdf_max_g, g_max] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2);
    pkf_value = pkf_3d(s1, s2,true,false);
    prob_g_max = pkf_value*pdf_max_g;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
    mu_vee = get_vee_vector(mu_g);
    s1_vee = zeros(6,1);
    pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
    prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
   
   % position error distribution of s2 (the case only consider the
    % position error but the real case has orientation error)
    mu = s2.tc;
    Sigma = zeros(3);
    Sigma(1,1) = 1.0000e-06;
    Sigma(2,2) = 4.0000e-06;
    Sigma(3,3) = 9.0000e-06;

    % Compare with the methods that only consider the position error
    prob_max_poly = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2);
    
     % For non-sphere objects, use their minkowskisum to get the volume
    m1 = s1.GetGradientsCanonical();
    mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
    x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

    %mink_shp = alphaShape(x_mink(1,21:380)', x_mink(2,21:380)', x_mink(3,21:380)');
    %mink_vol = volume(alphaShape);
    [mink_conv,mink_vol] = convhull(x_mink(1,:)', x_mink(2,:)', x_mink(3,:)');   
    
     % Compute prob_x_s1
     prob_x_s1 = mvnpdf(s1.tc, mu, Sigma) * mink_vol;
     
     % Compute exact_prob_x
    exact_prob_x = exact_contact_probability_pure_translation(mu, Sigma, s1, s2, 1e+04);

    result(1,i) = double(flag);
    result(2,i) = dist;
    result(3,i) = prob_g_max;
    result(4,i) = exact_prob_g;
    result(5,i) = prob_max_poly;
    result(6,i) = prob_x_s1;
    result(7,i) = exact_prob_x;
    
end

%% Visualize 
figure; hold on;
s1.PlotShape('b', 0.3);
s2.PlotShape('g', 0.3);
[pdf_max_g, g_max] = max_contact_probability_SE3(mu_g, Sigma_g, s1, s2);
 
s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

s3.tc = g_max(1:3,4);
s3.q = rotm2quat(g_max(1:3,1:3));
s3.PlotShape('r', 0.3);
pause(0.5)
n=100;
gi_rand_vee = mvnrnd(get_vee_vector(mu_g), Sigma_g, n);
for i = 1:n
    gi_rand_matrix = get_SE3_matrix(gi_rand_vee(i,1:6)');
    s3.tc = gi_rand_matrix(1:3,4);
    s3.q = rotm2quat(gi_rand_matrix(1:3,1:3));
    s3.PlotShape('r', 0.2);
end



