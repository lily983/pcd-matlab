 function [prob_max, g_max] = max_contact_probability_SO3_update_sigma(mu, Sigma, s1, s2, isplot)
% max_contact_probability Compute maximum probability density function value on the contact space
% between two convex bodies using closed-form parametric contact space and
% manifold optimization
%   
% Inputs
%   mu: Mean pose (center of s2)
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%Outputs
%   prob_max: The maximum pdf value 
%   g_max: 4by4 matrix of SE(3) group that has prob_max

if nargin < 5
    isplot = false;
end

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

% Create the problem structure.
so3 = rotationsfactory(3, 1);
problem.M = so3;
problem.cost = @(R) cost_function(R, mu, Sigma, s1, s3);

% Use auto-differentiation for gradient computations
% problem = manoptAD(problem);

% Solve the problem
options.maxiter = 100;
options.tolgradnorm = 1e-3;

[R_opt, f_opt, info] = trustregions(problem, [], options);

% Compute the max contact probability
prob_max = 1/( 8*pi^3 * sqrt(det(Sigma)) ) * exp(-1/2 * f_opt);

s3.q = rotm2quat(R_opt);
[~, mink_point] = dist_cfc_update_Sigma(s1, s3, Sigma(4:6,4:6));

p_opt = mink_point;

g_max = [R_opt, p_opt; 0, 0, 0, 1];

% Plots
if isplot
%     figure;
%     semilogy([info.iter], [info.gradnorm], '.-');
%     xlabel('Iteration number');
%     ylabel('Norm of the gradient of f');
    
    figure;
    semilogy([info.iter], [info.cost], '.-');
    xlabel('Iteration number');
    ylabel('Cost function value');
end
end

%% Objective function
function f = cost_function(R, mu, Sigma, s1, s3)
% Find closest points between s1 and s2
%s3.q = rotm2quat(R);

%{
[flag, dist, pt_cls] = distance_cfc(s1, s3);
p = pt_cls.mink;

if flag
    p = s3.tc;
end
%}

%{
opt =  optimoptions("fsolve","OptimalityTolerance",1e-15);
fun = @(lamda)norm((inv(Sigma(4:6,4:6)) + lamda*eye(3))\(Sigma(4:6,4:6)\mu(1:3,4) + lamda*s1.tc) - s1.tc) - s1.a(1) - s3.a(1);
lamda0 = 0.0;
lamda = fsolve(fun, lamda0, opt);

x_max =double((inv(Sigma(4:6,4:6)) + lamda*eye(3))\(Sigma(4:6,4:6)\mu(1:3,4)  + lamda*s1.tc));
%}

% Update s3 rotation, and deform Sigma to identity and find the new closed point on minksum
s3.q = rotm2quat(R);
[~, mink_point] = dist_cfc_update_Sigma(s1, s3, Sigma(4:6,4:6));

x_max = mink_point;

% Compute cost
g = [R, x_max; 0, 0, 0, 1];
x_diff = get_vee_vector(mu \ g);

f = x_diff' / Sigma * x_diff;

end