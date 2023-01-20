%% Test 1: pure traslation case. prob_x_max, prob_x_s1, prob_x_polygon, exact_prob_x
close all;clear;clc;
add_path()

N = [20,20];
flag  =1;

for i=1:50
    
    while flag
        s1 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
        zeros(3,1), [1, 0, 0, 0], N});

        s2 = SuperQuadrics({0.05*ones(1,3), [1,1], [0, 0]...
            rand(3,1), [1, 0, 0, 0], N});

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
mu=s2.tc;
Sigma = zeros(3);
Sigma(1,1) = 0.0001;
Sigma(2,2) = 0.0009;
Sigma(3,3) = 0.0025;
xi_rand = mvnrnd(mu, Sigma, n);
for i = 1:size(xi_rand, 1)
   
    s2.tc = xi_rand(i,1:3)';

    s2.PlotShape('r', 0.1);
end
