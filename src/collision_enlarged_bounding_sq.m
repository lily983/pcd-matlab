function flag = collision_enlarged_bounding_sq(s1, s2, mu, Sigma)
% This function computes the collision status using enlarged bounding
% convex mesh of s2 based on its pose error propability distribution
Patch_s2 = enlarged_bounding_superquadrics(s2, mu, Sigma);
patch_s2.vertices = Patch_s2.Vertices;
patch_s2.faces = Patch_s2.Faces;

s1_points = s1.GetPoints()';
patch_s1 = surf2patch(reshape(s1_points(:,1), [20,20]), reshape(s1_points(:,2), [20,20]), reshape(s1_points(:,3), [20,20]), 'triangles');

flag = GJK(patch_s1, patch_s2, 1000);

figure; hold on;
s1.PlotShape('b', 0.2)
s2.PlotShape('g', 0.2)
% bounding_s2.PlotShape('y',0.5)
% Patch_s2
patch(patch_s1, 'FaceAlpha', 0.3)
patch(patch_s2, 'FaceAlpha', 0.3)

end