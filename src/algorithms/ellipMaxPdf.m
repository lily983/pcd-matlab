function [prob,time]=ellipMaxPdf(e1, e2, xx, Sigmax)
%ellipMaxPdf: revised implementation of paper "Efficient Probabilistic
% Collision Detection for Non-Convex Shapes". The original paper is only 
% for two spherical objects, here we extend to ellipsoidal objects
%
%Inputs
%   e1, e2: Two ellipsoidal objects
%   xx: Mean of the relative position error x = x2-x1 
%   Sigmax: Covariance of the relative position error
%Outputs
%   prob: The probability approximation
%   time: computation time

% Start record algorithm running time
tic;
prob = 0;

% First do space transformation to make the bounding ellipsoid as a
% unit ball located
% Noted that the bounding ellipsoid is an upper bound for the Minkowski sum
% of two ellipsoids
[~, Sigmaf] = get_bounding_ellip(e1, e2);

% Do the same space transformation to the error distribution
new_xx = (sqrtm(Sigmaf)) \ xx;
new_Sigmax =  (sqrtm(Sigmaf)) \  Sigmax / sqrtm(Sigmaf);


% Inputs: new_xx (n×1), new_Sigmax (n×n, symmetric positive definite)
W = inv(new_Sigmax);

% Find the surface point on the boundary of the bounding ellipsoid, which has 
%the maximum pdf value. 
% Noted that here the bounding ellipsoid is transformed to be an sphere,
%so we need to normalize the point
cost_y = @(y) (y/norm(y) - new_xx).' * W * (y/norm(y) - new_xx);

% Initial value
y0 = new_xx;

options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off'); % no gradient needed

% Result
yopt = fminunc(cost_y, y0, options);

% project the result back onto the sphere
xopt = yopt / norm(yopt);              

pdf = mvnpdf(xopt.', new_xx.', new_Sigmax);

% The probability approximation is the pdf value times with the integration
% volume
prob = pdf * pi * 4/3;

time = toc;
end