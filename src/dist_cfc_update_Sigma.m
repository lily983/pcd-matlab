function [dist, surf_point] = dist_cfc_update_Sigma(s1, s2, Sigma)

m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s2,  quat2rotm(s1.q), quat2rotm(s2.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

% Transfome points so that the position error distribution are spheres
% centered at s2 center, N(x; s2, eye(3))
mink_points = (sqrtm(Sigma)\x_mink)';
shift = eye(3) * s2.tc - sqrtm(Sigma) \ s2.tc;

mink_points(:,1) = mink_points(:,1) + shift(1,1);
mink_points(:,2) = mink_points(:,2) + shift(2,1);
mink_points(:,3) = mink_points(:,3) + shift(3,1);

X_ = reshape(mink_points(:,1), [20,20]); 
Y_= reshape(mink_points(:,2),  [20,20]); 
Z_ = reshape(mink_points(:,3),  [20,20]); 

patch_mink = surf2patch(X_,Y_,Z_, 'triangles');

arg.Faces = patch_mink.faces;
arg.Vertices = patch_mink.vertices;
arg.QueryPoints = s2.tc';

[~, surf_point, ~, ~, ~, ~] = point2trimesh(arg);

% Transform points back 
surf_point = sqrtm(Sigma)*surf_point' + (eye(3) - sqrtm(Sigma))*s2.tc;

dist = norm(surf_point - s2.tc);
end