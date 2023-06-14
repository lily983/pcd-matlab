%% Test 1: SE(3) pose error case. prob_g_max, prob_g_s1, exact_prob_g, exact_prob_x, prob_x_poly, prob_x_s1
close all;clear;clc;
add_path()

N = [20,20];

s1 = SuperQuadrics({0.02*ones(1,3), [1.0, 0.8], [0, 0]...
zeros(3,1), [1, 0, 0, 0], N});

% col:40
for i=1:5
    flag = 1;
    while flag
        init_quat = [0.8678    0.3188    0.1610    0.3454] + rand(1,4)*0.01;
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
            [ 0.0145;0.0382;0.0366] + rand(3,1)*0.01,  init_quat, N}); % col:40
%         
%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   
    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result(1,i) = double(flag);
    result(2,i) = dist;
    result(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result(5,i) = prob_max_g_so3_cls;
    result(6,i) = exact_prob_g;
    result(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result(22:27,i) = get_vee_vector(g_max_so3_cls);
    result(29:34,i) = get_vee_vector(mu_g);
    result(36, i) = exact_prob_g_GJK;
    result(37,i) = prob_bounding_sq;
    
end

%col 36
for i=1:5
    flag = 1;
    while flag
        init_quat = rand(1,4)*0.01 + [0.3629    0.5736    0.7341    0.0213];
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
                     [0.0105; 0.0370; 0.0122] + rand(3,1)*0.01,  init_quat, N}); 

%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   

    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result2(1,i) = double(flag);
    result2(2,i) = dist;
    result2(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result2(5,i) = prob_max_g_so3_cls;
    result2(6,i) = exact_prob_g;
    result2(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result2(22:27,i) = get_vee_vector(g_max_so3_cls);
    result2(29:34,i) = get_vee_vector(mu_g);
    result2(36, i) = exact_prob_g_GJK;
    result2(37,i) = prob_bounding_sq;
    
end

%col 35
for i=1:5
    flag = 1;
    while flag
        init_quat = rand(1,4)*0.01 + [0.7598    0.4160    0.4158    0.2770];
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
            [0.0012; 0.0397; 0.0286] + rand(3,1)*0.01, init_quat, N});
        
%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   

    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result3(1,i) = double(flag);
    result3(2,i) = dist;
    result3(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result3(5,i) = prob_max_g_so3_cls;
    result3(6,i) = exact_prob_g;
    result3(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result3(22:27,i) = get_vee_vector(g_max_so3_cls);
    result3(29:34,i) = get_vee_vector(mu_g);
    result3(36, i) = exact_prob_g_GJK;
    result3(37,i) = prob_bounding_sq;
    
end

%col 22
for i=1:5
    flag = 1;
    while flag
        init_quat = rand(1,4)*0.01 + [0.3004    0.5080    0.6213    0.5155];
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
%         s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
%             init_posi, init_quat, N});
      
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
            [0.0338; 0.0113; 0.0143] + rand(3,1)*0.01, init_quat, N});
        
%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   

    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result4(1,i) = double(flag);
    result4(2,i) = dist;
    result4(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result4(5,i) = prob_max_g_so3_cls;
    result4(6,i) = exact_prob_g;
    result4(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result4(22:27,i) = get_vee_vector(g_max_so3_cls);
    result4(29:34,i) = get_vee_vector(mu_g);
    result4(36, i) = exact_prob_g_GJK;
    result4(37,i) = prob_bounding_sq;
    
end

%col 16
for i=1:5
    flag = 1;
    while flag
        init_quat = rand(1,4)*0.01+[0.5712    0.6756    0.2715    0.3789];
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
%         s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
%             init_posi, init_quat, N});
      
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
            [0.0011; 0.0245; 0.0341] + rand(3,1)*0.01, init_quat, N});
        
%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   

    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result5(1,i) = double(flag);
    result5(2,i) = dist;
    result5(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result5(5,i) = prob_max_g_so3_cls;
    result5(6,i) = exact_prob_g;
    result5(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result5(22:27,i) = get_vee_vector(g_max_so3_cls);
    result5(29:34,i) = get_vee_vector(mu_g);
    result5(36, i) = exact_prob_g_GJK;
    result5(37,i) = prob_bounding_sq;
    
end

%col 14
for i=1:5
    flag = 1;
    while flag
        init_quat = rand(1,4)*0.01 + [0.6911    0.3981    0.6009    0.0532];
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
%         s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
%             init_posi, init_quat, N});
      
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
            [0.0137; 0.0464; 0.0247] + rand(3,1)*0.01, init_quat, N});
        
%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   

    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result6(1,i) = double(flag);
    result6(2,i) = dist;
    result6(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result6(5,i) = prob_max_g_so3_cls;
    result6(6,i) = exact_prob_g;
    result6(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result6(22:27,i) = get_vee_vector(g_max_so3_cls);
    result6(29:34,i) = get_vee_vector(mu_g);
    result6(36, i) = exact_prob_g_GJK;
    result6(37,i) = prob_bounding_sq;
    
end

%col 10
for i=1:5
    flag = 1;
    while flag
        init_quat = rand(1,4)*0.01 + [0.7305    0.6675    0.1443    0.0020];
        init_quat = init_quat/norm(init_quat);
%         init_posi = 0.035*rand(3,1);

%         [0.024058038057378;0.026986650282266;0.026057953781762]
%         s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
%             init_posi, init_quat, N});
      
        s2 = SuperQuadrics({[0.01,0.02,0.01], [1,0.8], [0, 0]...
            [0.0340; 0.0220; 0.0382] + rand(3,1)*0.01, init_quat,  N});
        
%         figure; hold on;
%         s1.PlotShape('b', 0.3);
%         s2.PlotShape('g', 0.3);
%         pause(1);
        
        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    % pose error distribution of s2
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g = zeros(6);
    Sigma_g(1,1) =1.0000e-2;
    Sigma_g(2,2) =1.0000e-2;
    Sigma_g(3,3) =1.0000e-2;
    Sigma_g(4,4) =1.0000e-5;
    Sigma_g(5,5) =1.0000e-5;
    Sigma_g(6,6) =1.0000e-5;
   

    pkf_value = pkf_3d(s1, s2,true,false);
    
    try
    % Compute prob_g_max on manifold SO3*S2
        [pdf_max_g_so3_s2, g_max_so3_s2] = max_contact_probability_SO3_S2(mu_g, Sigma_g, s1, s2,false);
    catch 
        pdf_max_g_so3_s2 = NaN;
    end
    prob_max_g_so3_s2 = pkf_value*pdf_max_g_so3_s2;
    
    % Compute prob_g_max on manifold SO3 with x_max
%     [pdf_max_g_so3_xmax, g_max_so3_xmax] = max_contact_probability_SO3_lamda_search(mu_g, Sigma_g, s1, s2,false);
%     prob_max_g_so3_xmax = pkf_value*pdf_max_g_so3_xmax;
    
    % Compute prob_g_max on manifold SO3 with closed point
    try
        [pdf_max_g_so3_cls, g_max_so3_cls] = max_contact_probability_SO3_update_sigma(mu_g, Sigma_g, s1, s2,false);
    catch
        pdf_max_g_so3_cls = NaN;
    end
    prob_max_g_so3_cls = pkf_value*pdf_max_g_so3_cls;
   
    % Compute prob_g_s1, get pdf_g value of error distribution of s2 at the pose of s1 
%     mu_vee = get_vee_vector(mu_g);
%     s1_vee = zeros(6,1);
%     pdf_g_s1 = mvnpdf(s1_vee', mu_vee', Sigma_g);
%     prob_g_s1 = pkf_value*pdf_g_s1;
    
    % Compute exact prob_g 
    try
         exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g = NaN;
    end

     % Compute exact probability by using patch based GJK algorithm
    try 
         exact_prob_g_GJK = exact_contact_probability_collision_mesh(mu_g, Sigma_g, s1, s2, 1e+04);
    catch
        exact_prob_g_GJK = NaN;
    end

     % Compute the collision status using enlarged bounding sq
     try
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
    end

    result7(1,i) = double(flag);
    result7(2,i) = dist;
    result7(3,i) = prob_max_g_so3_s2;
%     result(4,i) = prob_max_g_so3_xmax;
    result7(5,i) = prob_max_g_so3_cls;
    result7(6,i) = exact_prob_g;
    result7(8:13,i) = get_vee_vector(g_max_so3_s2);
%     result(15:20,i) = get_vee_vector(g_max_so3_xmax);
    result7(22:27,i) = get_vee_vector(g_max_so3_cls);
    result7(29:34,i) = get_vee_vector(mu_g);
    result7(36, i) = exact_prob_g_GJK;
    result7(37,i) = prob_bounding_sq;
    
end

%% Visualize 

init_quat = rand(1,4);
init_quat = init_quat/norm(init_quat);
init_posi = 0.045*rand(3,1);

s2 = SuperQuadrics({[0.01,0.02,0.03], [0.4,0.8], [0, 0]...
    init_posi, init_quat, N});
%%
figure; hold on;
s1.PlotShape('b', 0.3);
s2.PlotShape('r', 0.1);
%%
m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

surf(X, Y, Z,...
 'EdgeColor', 'k', 'EdgeAlpha', 0.2,...
 'FaceAlpha', 0.1, 'FaceColor', 'b'); %plot contact space

[flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'fixed-point');

scatter3(pt_cls.mink(1),pt_cls.mink(2),pt_cls.mink(3),'r','+');


%%
s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});
mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
 Sigma_g = zeros(6);
% Sigma_g(1,1) =1.0000e-06;
% Sigma_g(2,2) =9.0000e-06;
% Sigma_g(3,3) =8.0000e-06;
% Sigma_g(4,4) =1.0000e-06;
% Sigma_g(5,5) =4.0000e-06;
% Sigma_g(6,6) =9.0000e-06;
% 
%     Sigma_g(1,1) =1.0000e-6;
%     Sigma_g(2,2) =1.0000e-6;
%     Sigma_g(3,3) =1.0000e-6;
%     Sigma_g(4,4) =1.0000e-6;
%     Sigma_g(5,5) =1.0000e-6;
%     Sigma_g(6,6) =1.0000e-6;
% 
    Sigma_g(1,1) =1.0000e-3;
    Sigma_g(2,2) =1.0000e-3;
    Sigma_g(3,3) =1.0000e-3;
    Sigma_g(4,4) =1.0000e-6;
    Sigma_g(5,5) =1.0000e-6;
    Sigma_g(6,6) =1.0000e-6;

n=200;
gi_rand_vee = mvnrnd(get_vee_vector(mu_g), Sigma_g, n);
for i = 1:n
    gi_rand_matrix = get_SE3_matrix(gi_rand_vee(i,1:6)');
    s3.tc = gi_rand_matrix(1:3,4);
    s3.q = rotm2quat(gi_rand_matrix(1:3,1:3));
    s3.PlotShape('r', 0.05);

end

%% Store results
pathname = '~/prob-collision-matlab/test/result';
test = "se3_R_10_2_t_10_5_cfc_CS_col10" ;
filename = test + ".mat";
save(fullfile(pathname, filename), 'result7')
