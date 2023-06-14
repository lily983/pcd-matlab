function visualize_3D_case(result, col, Sigma,flag)

%flag=1: ellip, flag=0, sphere
if flag==0
    s1 = SuperQuadrics({result(11,col)*ones(1,3), [1,1], [0, 0]...
    result(12:14,col), result(15:18,col)', [20, 20]});

    s2 = SuperQuadrics({result(20,col)*ones(1,3),  [1,1], [0, 0]...
       result(21:23,col), result(24:27,col)',[20,20]});
elseif flag==1
    s1 = SuperQuadrics({result(11:13,col)', [1,1], [0, 0]...
    result(14:16,col), result(17:20,col)', [20, 20]});

    s2 = SuperQuadrics({result(22:24,col)',  [1,1], [0, 0]...
       result(25:27,col), result(28:31,col)',[20,20]});
end

visualize_bounding_ellip(s1, s2);

visualize_position_error(s1, s2, Sigma)

% % Get points that gmm(x)<iota(x)
% test_gmm_3D(s1, s2) % test if GMM is upper bound

view(0,75)
lightangle(-45,30)
% h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.1;
h.DiffuseStrength = 0.3;
h.SpecularStrength = 0.3;
h.SpecularExponent = 15;
% h.BackFaceLighting = 'unlit';

end