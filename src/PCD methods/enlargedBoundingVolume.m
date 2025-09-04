function [prob, time] = enlargedBoundingVolume(s1, s2, Sigma1, Sigma2, confidenceLevel, isPlot)
% enlargedBoundingVolume: Enlarge the object s1 and s2 by the 
% uncertainty ellipsoid constructed by Sigma1 and Sigma2 and
% confidenceLevel (1-epsilon, e.g. 99%confidenceLevel->epsilon=0.01)
% This function is the implementation of paper "Fast Certification of
% Collision Probability Bounds with Uncertain Convex Obstacles"
% The original paper only consider one object has position uncertainty, in
% this function, we use Eq. 7 in the paper to construct two enlarged
% objects based on their covariance matrix and confidenceLevel. Then, use
% GJK to check if they collide.
    % Inputs:
        % s1, s2: sphere or ellipsoid or superquadraics
        % Sigma1: covariance matrix of position error distribution of s1
        % Sigma2: covariance matrix of position error distribution of s2
        % confidenceLevel: 1-epsilon
    % Outputs:
        % prob: PCD value
        % t: computation time
if nargin == 5
    isPlot = false;
end

prob = 0;
tic;

% check if Sigma is zeros
if ~all(diag(Sigma1))
    s1_enlarged_points = s1.GetPoints();
else
    % position error ellipsoid
    scale = chi2inv(confidenceLevel, 3);
    % here error_ellip is assumed to be diagonal matrix defined in the
    % world space, so it has identity rotation matrix
    error_ellip = SuperQuadrics({(scale.*diag(Sigma1)).^0.5, [1,1], [0,0], [0;0;0], [1,0,0,0], s1.N});
    s1_enlarged = MinkSumClosedForm(s1, error_ellip, quat2rotm(s1.q), quat2rotm(error_ellip.q));
    m1 = s1.GetGradientsCanonical();
    s1_enlarged_points = s1_enlarged.GetMinkSumFromGradient(m1) + s1.tc;
end

if ~all(diag(Sigma2))
    s2_enlarged_points = s2.GetPoints();
else
    % position error ellipsoid
    scale = chi2inv(confidenceLevel, 3);
    error_ellip2 = SuperQuadrics({(scale.*diag(Sigma2)).^0.5, [1,1], [0,0], [0;0;0], [1,0,0,0], s2.N});
    s2_enlarged = MinkSumClosedForm(s2, error_ellip2, quat2rotm(s2.q), quat2rotm(error_ellip2.q));
    m2 = s2.GetGradientsCanonical();
    s2_enlarged_points = s2_enlarged.GetMinkSumFromGradient(m2) + s2.tc;
end

SN1 = s1.N;
s1_enlarged_patch = surf2patch(reshape(s1_enlarged_points(1,:),...
    SN1),reshape(s1_enlarged_points(2,:),  SN1),reshape(s1_enlarged_points(3,:),  SN1), 'triangles');

SN2 = s2.N;
s2_enlarged_patch = surf2patch(reshape(s2_enlarged_points(1,:),...
    SN2),reshape(s2_enlarged_points(2,:),  SN2),reshape(s2_enlarged_points(3,:),  SN2), 'triangles');

prob = GJK(s1_enlarged_patch, s2_enlarged_patch, 1000);
time = toc;

if isPlot
%     visualize_bounding_ellip(s1, s2);
    % 
    s1_enlarged_patch.FaceColor = 'y';
    s1_enlarged_patch.FaceAlpha = 0.5;
    patch(s1_enlarged_patch);
    % 
    s2_enlarged_patch.FaceColor = 'm';
    s2_enlarged_patch.FaceAlpha = 0.5;
    patch(s2_enlarged_patch);
    % % 
%     visualize_bounding_ellip(s1, s2);
%     error_ellip2.tc = s2.tc;
%     error_ellip.PlotShape('y', 0.5);
%     error_ellip2.PlotShape('m', 0.5);
end
end