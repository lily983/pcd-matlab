function flag = collision_mesh(s1, s2)
% This function compute the collision states of two meshed objects using
% GJK algorithm

s1_points = s1.GetPoints()';

s2_points = s2.GetPoints()';

dimension = size(s1.a, 2);

if dimension == 2
    flag = collisionCheck2D(s1_points, s2_points);
elseif dimension == 3
    patch_s1 = surf2patch(reshape(s1_points(:,1), s1.N), reshape(s1_points(:,2), s1.N), reshape(s1_points(:,3), s1.N), 'triangles');
    patch_s2 = surf2patch(reshape(s2_points(:,1), s2.N), reshape(s2_points(:,2), s2.N), reshape(s2_points(:,3), s2.N), 'triangles');
    flag = GJK(patch_s1, patch_s2, 1000);
end

% figure; hold on;
%  patch(patch_s1);
% patch(patch_s2);
% [dist,pts,G,H] = GJK_dist(patch_s1, patch_s2)

end