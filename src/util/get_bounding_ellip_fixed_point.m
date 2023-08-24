function [mf, Sigmaf] = get_bounding_ellip_fixed_point(s1, s2)
% This function gets the bounding ellip for the Minkowski sum of s1 and s2
% Get beta by fixed-point iteration

% Get defining matrix of ellipsoid
Q1 = quat2rotm(s1.q) * diag(s1.a).^2 * quat2rotm(s1.q)';
Q2 = quat2rotm(s2.q) * diag(s2.a).^2 * quat2rotm(s2.q)';

% Initial value of optimization
beta0 = sqrt(trace(Q1)/trace(Q2));

% Optimization
beta_opt = fixed_point(beta0, Q1, Q2);

Sigmaf = (1 + 1.0/beta_opt) * Q1 + (1 + beta_opt) * Q2;
mf = s1.tc - s2.tc;

end

% Fixed-point method
function beta_new = fixed_point(beta0, Q1, Q2)
% Tolerance
tol = 1e-10;

% Lamda
R = Q1 \ Q2;
lamda = eig(R);

beta_initial = beta0;

for i=1:500
    beta_new = fixed_point_iteration(beta_initial, lamda);
    
    if abs(beta_new - beta0) < tol
      return;
    end
    beta_initial = beta_new;
end

end

% Fixed-point iteration
function beta1 = fixed_point_iteration(beta0, lamda)

beta1 = sqrt( ( 1/(1+beta0*lamda(1)) + 1/(1+beta0*lamda(2)) + 1/(1+beta0*lamda(3)) ) ...
    /  ( lamda(1)/(1+beta0*lamda(1)) + lamda(2)/(1+beta0*lamda(2)) + lamda(3)/(1+beta0*lamda(3)) ) );
end