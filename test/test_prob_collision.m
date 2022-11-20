close all; clear; clc;
add_path()

N = [20,20];
s1 = SuperQuadrics({0.5+5*rand(1,3), 0.01+1.98*rand(1,2), [0,0],...
    zeros(3,1), [1,0,0,0], N});
s2 = SuperQuadrics({0.5+5*rand(1,3), 0.01+1.98*rand(1,2), [0,0],...
    5*(2*rand(3,1)-1), rand(1,4), N});

Sigma = 0.1*eye(6);
Sigma(4:6,4:6) = 0.5*eye(3);

%% Solve max contact probability
mu = [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];

[prob_max, g_max] = max_contact_probability(mu, Sigma, s1, s2, true);
pkf_value = pkf_3d(s1, s2, true, false);
prob_max = pkf_value * prob_max;

%% Plots
figure; hold on; axis equal;
lightangle(gca,45,30);
lighting gouraud;

s1.PlotShape('b', 0.7);
s2.PlotShape('g', 0.7);

trplot(g_max)

% Sample noisy s2 poses
n = 100;

mu_log_skew = logm(mu);
mu_log = [vex(mu_log_skew(1:3,1:3)); mu_log_skew(1:3,4)];

xi_rand = mvnrnd(mu_log', Sigma, n);
for i = 1:size(xi_rand, 1)
    xi_skew = [skew(xi_rand(i,1:3)), xi_rand(i,4:6)'; zeros(1,4)];
    g_rand = expm(xi_skew);
    
    s2.tc = g_rand(1:3,4);
    s2.q = rotm2quat(g_rand(1:3,1:3));
    s2.PlotShape('r', 0.1);
end