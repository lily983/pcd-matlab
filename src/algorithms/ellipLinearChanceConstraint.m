function [prob, t, a, x_mink] = ellipLinearChanceConstraint(s1, s2, xx, Sigma)
% This function reproduce methods using linear chance constraint (LCC)
% to get PCD value for Gaussian distributed variable. The code is inspired
% based on papers "Chance-Constrained Collision Avoidance for MAVs in
% Dynamic Environments" and "Tight Collision Probability for UAV Motion
% Planning in Uncertain Environment" 
% Inputs:
% s1, s2: sphere or ellipsoid or superquadraics
% xx: Mean of relative position error x = x2-x1
% Sigma: covariance matrix of position error distribution
% Outputs:
% prob: PCD value
% t: computation time

tic;
prob = 0;

[~, Sigmaf] = get_bounding_ellip(s1, s2);
xx = (s2.tc - s1.tc);
a = (Sigmaf^0.5 \ xx) ./ norm(Sigmaf^0.5 \ xx);
prob = 1/2 + 1/2 * erf( (1-a' / Sigmaf^0.5 * xx) / sqrt(2*a' / Sigmaf^0.5 * Sigma / Sigmaf^0.5 *a));
t = toc;
x_mink = xx ./ norm(Sigmaf^0.5 \ xx);
a = (Sigmaf \ xx) ./ norm(Sigmaf \ xx);

end