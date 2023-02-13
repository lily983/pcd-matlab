function   [pdf_max_x, x_max] = max_contact_probability_pure_translation_superquadrics(mu, Sigma, s1, s2)
% max_probability_pure_translation: Center position of s2 is a Guassian random variable, 
% the function calculate the maximum probability that it will collide with s1 
%
%Inputs
%   mu: Mean pose (center of s2)
%   Sigma: Covariance of the probability
%   s1, s2: SuperQuadrics
%Outputs
%   x_max: Center position of s2 that has the maximum probability of
%   colliding with s1
%   pdf_max_x: The probability of x_max

% Function file shares the same pointer to s2 in the main file
s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

if collision_cfc(s1,s3)
    x_max = s3.tc;
    pdf_max_x  = mvnpdf(x_max, mu, Sigma);
    return
end

mink = MinkSumClosedForm(s1, s3, quat2rotm(s1.q), quat2rotm(s3.q));

R2 = euclideanfactory(2);

problem.M = R2;
problem.cost = @(phi)cost_function(phi, s1, s3, mink, Sigma);

%{
% Initial point of optimization
R1 = quat2rotm(s1.q);
s3_tc_in_s1 = R1' * (s3.tc-s1.tc);
phi0 = [atan2( s3_tc_in_s1(3), norm(s3_tc_in_s1(1:2)) ),...
            atan2( s3_tc_in_s1(2), s3_tc_in_s1(1) )];
%}

% Solve the problem
options.maxiter = 100;
options.tolgradnorm = 1e-2;

[phi_opt, f_opt, ~] = trustregions(problem, [], options);

m1_opt = s1.GetGradientsFromSpherical(phi_opt');
x_max = mink.GetMinkSumFromGradient(m1_opt);

pdf_max_x = 1/( 8*pi^3 * sqrt(det(Sigma)) ) * exp(-1/2 * f_opt);

end

function f = cost_function(phi, s1, s3, mink, Sigma)

m1 = s1.GetGradientsFromSpherical(phi');
mink_point = mink.GetMinkSumFromGradient(m1);

f = (mink_point - s3.tc)' / Sigma * (mink_point - s3.tc);

end