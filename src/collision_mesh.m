function flag = collision_mesh(s1, s2)
% This function compute the collision states of two meshed objects using
% GJK algorithm

s1_points = s1.GetPoints()';

s2_points = s2.GetPoints()';

patch_s1 = surf2patch(reshape(s1_points(:,1), [20,20]), reshape(s1_points(:,2), [20,20]), reshape(s1_points(:,3), [20,20]), 'triangles');
patch_s2 = surf2patch(reshape(s2_points(:,1), [20,20]), reshape(s2_points(:,2), [20,20]), reshape(s2_points(:,3), [20,20]), 'triangles');

% figure; hold on;
% patch(patch_s1);
% patch(patch_s2);

% [dist,pts,G,H] = GJK_dist(patch_s1, patch_s2);
flag = GJK(patch_s1, patch_s2, 1000);

end