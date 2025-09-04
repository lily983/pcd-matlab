function [prob_max, g_max] = max_contact_probability_SE3(mu, Sigma, s1,s2)

% Function file shares the same pointer to s2 in the main file
s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

% Check if s1 and s3 collide at s3 mean pose, if collide, then the maximum
% probability is at the mean pose
[flag, ~, ~] = collision_cfc(s1, s3);

if flag
    g_max = mu;
    mu_vee = get_vee_vector(mu);
    g_max_vee = mu_vee;
    prob_max = mvnpdf(g_max_vee', mu_vee', Sigma);
    return
end

SE3 = specialeuclideanfactory(3,1);
problem.M = SE3;
problem.cost = @(g)cost_function(g, mu,Sigma);

% Solve the problem
options.maxiter = 100;
options.tolgradnorm = 1e-3;

[g_opt, f_opt, info] = trustregions(problem, [], options);

g_max = [g_opt.R, g_opt.t; 0,0,0,1];
prob_max = 1/( 8*pi^3 * sqrt(det(Sigma)) ) * exp(-1/2 * f_opt);

end

function f = cost_function(g, mu, Sigma)

g_matrix = [g.R, g.t; 0,0,0,1];
x_diff = get_vee_vector(mu \ g_matrix);
f = x_diff' / Sigma * x_diff;

end