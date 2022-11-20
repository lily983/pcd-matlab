function [prob_max, g_max] = max_contact_probability(mu, Sigma, s1, s2, isplot)
% max_contact_probability Compute maximum probability on the contact space
% between two convex bodies using closed-form parametric contact space and
% manifold optimization
%   
% Inputs
%   mu: Mean pose (center of s2)
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies

if nargin < 5
    isplot = false;
end

% Create the problem structure.
so3 = rotationsfactory(3, 1);
problem.M = so3;
problem.cost = @(R) cost_function(R, mu, Sigma, s1, s2);

% Use auto-differentiation for gradient computations
% problem = manoptAD(problem);

% Solve the problem
options.maxiter = 100;
options.tolgradnorm = 1e-3;

[R_opt, f_opt, info] = trustregions(problem, [], options);

% Compute the max contact probability
prob_max = 1/( 8*pi^3 * sqrt(det(Sigma)) ) * exp(-1/2 * f_opt);

% Get g_max
s2.q = rotm2quat(R_opt);
[~, ~, pt_cls] = distance_cfc(s1, s2);
p_opt = pt_cls.mink;

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
function f = cost_function(R, mu, Sigma, s1, s2)
% Find closest points between s1 and s2
s2.q = rotm2quat(R);
[flag, dist, pt_cls] = distance_cfc(s1, s2);
p = pt_cls.mink;

if flag
    p = s2.tc;
end

% Compute cost
g = [R, p; 0, 0, 0, 1];
x_diff_skew = logm(mu \ g);
x_diff = [vex(x_diff_skew(1:3,1:3)); x_diff_skew(1:3,4)];

f = x_diff' / Sigma * x_diff;

end