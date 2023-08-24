function visualize_3D_case(result1, result2, col, Sigma)

s1 = SuperQuadrics({result1(1:3,col)', [1,1], [0, 0]...
result1(8:10,col), result1(4:7,col)', [20, 20]});

s2 = SuperQuadrics({result2(1:3,col)',  [1,1], [0, 0]...
   result2(8:10,col), result2(4:7,col)',[20,20]});


visualize_bounding_ellip(s1, s2);

visualize_position_error(s1, s2, Sigma)

% figure; hold on; axis equal; axis off;
% view(0,75)
% lightangle(-45,30)
% % h.FaceLighting = 'gouraud';
% h.AmbientStrength = 0.1;
% h.DiffuseStrength = 0.3;
% h.SpecularStrength = 0.3;
% h.SpecularExponent = 15;
% h.BackFaceLighting = 'unlit';

end