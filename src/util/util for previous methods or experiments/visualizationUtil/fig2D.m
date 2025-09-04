close all;clear;clc;

% figure; hold on; axis equal; axis off;
figure; hold on; 

s1 = SuperEllipse([4, 4, 1, 0, 1, 6, 0.3, 25]);
s2 = SuperEllipse([5,5, 1, 0, 12, 5, 1.4, 25]);

% s1.PlotShape(hex2rgb('4FAAD1'), 0.7);
% s2.PlotShape(hex2rgb('45AC59'), 1);

Sigma = eye(2, 2) * 8*1e-02;
% visualize_position_error(s1, s2, Sigma)

visualize_bounding_ellip(s1, s2, Sigma)
% 
 [mf, Sigmaf] = get_bounding_ellip(s1, s2);
%  
 visualize_gmm_2D(mf, Sigmaf, 10, 3)

%%
for i=1:1
    flag  =1;
    while flag
        s1 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
        zeros(3,1), [1, 0, 0, 0], N});

        s2 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
           [0.104196729695369;0.038003107777440;0.006990404519385], [1, 0, 0, 0], N});

        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end

    mu=s2.tc;
    
    Sigma = zeros(3);
    Sigma(1,1) = 1.0000e-6;
    Sigma(2,2) = 4.0000e-6;
    Sigma(3,3) = 9.0000e-6;
    
   [pdf_max_x, x_max] = max_contact_probability_pure_translation(mu, Sigma, s1, s2);
 
    v =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
    prob_max_x = v * pdf_max_x;
    
    prob_x_s1 = mvnpdf(s1.tc, mu, Sigma) * v;

    exact_prob_x = exact_contact_probability_pure_translation(mu, Sigma, s1, s2, 1e+06);
  
    prob_max_poly = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2);
   
    result(1,i) = double(flag);
    result(2,i) = dist;
    result(3,i) = pdf_max_x;
    result(4,i) = prob_max_x;
    result(5,i) = prob_x_s1;
    result(6,i) = prob_max_poly;
    result(7,i) = exact_prob_x;
    result(9:11,i) = s2.tc;
    result(13:15,i) = x_max;

end    
%% Test 2. Pure translation error. Compare when sigma has different diagonal values and the same diagonal value
% When sigma has the same diagonal value: x_max is the same point as
% pt_cls.mink
% When sigma has the different diagonal value: x_max and pt_cls.mink are
% not the same point. PDF value at x_max is much larger than pt_cls.mink
N=[20,20]
s1 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
        zeros(3,1), [1, 0, 0, 0], N});

s2 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
    0.15*rand(3,1), [1, 0, 0, 0], N});

[flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'fixed-point');

mu=s2.tc;

Sigma_same = 3.0000e-05*eye(3);

Sigma_different = zeros(3);
Sigma_different(1,1) = 1.0000e-05;
Sigma_different(2,2) = 4.0000e-05;
Sigma_different(3,3) = 9.0000e-05;

m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

figure; hold on;
surf(X, Y, Z,...
 'EdgeColor', 'k', 'EdgeAlpha', 0.2,...
 'FaceAlpha', 0.1, 'FaceColor', 'b'); %plot contact space

 s1.PlotShape('b', 0.3);
 s2.PlotShape('g', 0.3);

% x_max when sigma has the same diagonal value
[pdf_max_x_ss, x_max_ss] = max_contact_probability_pure_translation(mu, Sigma_same, s1, s2)
scatter3(x_max_ss(1), x_max_ss(2), x_max_ss(3),	200, 'r');

% x_max when sigma has the different diagonal value
[pdf_max_x_sd, x_max_sd] = max_contact_probability_pure_translation(mu, Sigma_different, s1, s2)
scatter3(x_max_sd(1), x_max_sd(2), x_max_sd(3),200,"MarkerEdgeColor","#9932CC");

v =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
prob_max_x_sd = v*pdf_max_x_sd;

exact_prob_x_sd = exact_contact_probability_pure_translation(mu, Sigma_different, s1, s2, 1e+05);

scatter3(pt_cls.mink(1), pt_cls.mink(2), pt_cls.mink(3), 80,'b');
%% Store results
pathname = '~/prob-collision-matlab/test/result';
test = "compare_max_x_and_max_poly" + "_exclude_collide";
filename = test + ".mat";
fig = test +".fig";
save(fullfile(pathname, filename), 'result')
savefig(fullfile(pathname, fig))


%% Show the covariance effect on s2
n = 100;

xi_rand = mvnrnd(s2.tc, Sigma, n);

m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

figure; hold on;
s2.PlotShape('r', 0.1);
s1.PlotShape('g', 0.3);

surf(X, Y, Z,...
 'EdgeColor', 'k', 'EdgeAlpha', 0.1,...
 'FaceAlpha', 0.1, 'FaceColor', 'b'); %plot contact space

for i = 1:size(xi_rand, 1)
   
    s2.tc = xi_rand(i,1:3)';

    s2.PlotShape('r', 0.03);
end

x_max = result(13:15,3);

%scatter3(x_max(1), x_max(2), x_max(3),200,'g*');

%scatter3(s1.tc(1), s1.tc(2), s1.tc(3), 200, 'r');

[~, ~, pt_cls, ~] = collision_cfc(s1,s2,'fixed-point');

%scatter3(pt_cls.mink(1), pt_cls.mink(2), pt_cls.mink(3), 200,'b+');
