function [prob_max, g_max] = max_contact_probability_SO3_S2(mu, Sigma, s1,s2, isplot)

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

SO3_S2= productmanifold(struct('R', rotationsfactory(3), 'phi', spherefactory(2)));

problem.M = SO3_S2;
problem.cost = @(R_phi)cost_function(R_phi, s1, s3, mu,Sigma);


% Initial point of optimization
R1 = quat2rotm(s1.q);
R2 = quat2rotm(s2.q);
s2_tc_in_s1 = R1' * (s2.tc-s1.tc);
psi0 = [atan2( s2_tc_in_s1(3), norm(s2_tc_in_s1(1:2)) ),...
            atan2( s2_tc_in_s1(2), s2_tc_in_s1(1) )];
R_phi0.R = R2;
R_phi0.phi = psi0';


% Solve the problem
options.maxiter = 100;
options.tolgradnorm = 1e-4;

[R_phi_final, f_opt, info] = trustregions(problem, [], options);

mink_final = MinkSumClosedForm(s1, s3, quat2rotm(s1.q), R_phi_final.R);

% Find the final points on the minkowski sum
m1_final = s1.GetGradientsFromSpherical(R_phi_final.phi');
mink_poin_final = mink_final.GetMinkSumFromGradient(m1_final);

g_max = [R_phi_final.R, mink_poin_final; 0,0,0,1];
prob_max = 1/( 8*pi^3 * sqrt(det(Sigma)) ) * exp(-1/2 * f_opt);
% Plot

if nargin < 5
    isplot = false;
end

if isplot
    figure;
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the gradient of f');
    
    figure;
    semilogy([info.iter], [info.cost], '.-');
    xlabel('Iteration number');
    ylabel('Cost function value');
end

end

function f = cost_function(R_phi, s1, s3, mu, Sigma)

% Get minksum surface point from spherical parametrization
mink = MinkSumClosedForm(s1, s3, quat2rotm(s1.q), R_phi.R);

m1 = s1.GetGradientsFromSpherical(R_phi.phi');
mink_point = mink.GetMinkSumFromGradient(m1);

g_matrix = [R_phi.R, mink_point; 0,0,0,1];
x_diff = get_vee_vector(mu \ g_matrix);
f = x_diff' / Sigma * x_diff;

end