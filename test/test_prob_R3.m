%% Test 1: pure traslation case. prob_x_max, prob_x_s1, prob_x_polygon, exact_prob_x
close all;clear;clc;
add_path()

N = [20,20];

for i=1:50
    flag  =1;
    while flag
        s1 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
        zeros(3,1), [1, 0, 0, 0], N});

        s2 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
            0.2*rand(3,1), [1, 0, 0, 0], N});

        [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    end
    
    mu=s2.tc;
    
    Sigma = zeros(3);
    Sigma(1,1) = 1.0000e-06;
    Sigma(2,2) = 9.0000e-06;
    Sigma(3,3) = 2.5000e-05;
    
   [pdf_max_x, x_max] = max_probability_pure_translation(mu, Sigma, s1, s2);
 
    v =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
    prob_max_x = v * pdf_max_x;
    
    prob_x_s1 = mvnpdf(s1.tc, mu, Sigma) * v;

    exact_prob_x = exact_contact_probability_pure_translation(mu, Sigma, s1, s2, 1e+05);
  
    prob_max_poly = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2);
   
    result(1,i) = double(flag);
    result(2,i) = dist;
    result(3,i) = pdf_max_x;
    result(4,i) = prob_max_x;
    result(5,i) = prob_x_s1;
    result(6,i) = prob_max_poly;
    result(7,i) = exact_prob_x;
    
     s1.PlotShape('b', 0.5);
     s2.PlotShape('g', 0.5);

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
[pdf_max_x_ss, x_max_ss] = max_probability_pure_translation(mu, Sigma_same, s1, s2)
scatter3(x_max_ss(1), x_max_ss(2), x_max_ss(3),	200, 'r');

% x_max when sigma has the different diagonal value
[pdf_max_x_sd, x_max_sd] = max_probability_pure_translation(mu, Sigma_different, s1, s2)
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

%% Plot results
figure; hold on;

plot(result(7,:),'b')
plot(result(6,:),'g')

%% Box plot to compare upper bound and exact probability
figure;
difference  = result(4,:)' - result(3,:)';
boxplot(difference)

%% Show the covariance effect on s2
n = 100;
Sigma_different = zeros(3);
Sigma_different(1,1) = 1.0000e-05;
Sigma_different(2,2) = 4.0000e-05;
Sigma_different(3,3) = 9.0000e-05;
xi_rand = mvnrnd(mu, Sigma_different, n);
for i = 1:size(xi_rand, 1)
   
    s2.tc = xi_rand(i,1:3)';

    s2.PlotShape('r', 0.1);
end

