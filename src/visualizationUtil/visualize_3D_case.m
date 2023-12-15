function visualize_3D_case(result1, result2, col, Sigma1, Sigma2)

s1 = SuperQuadrics({result1(1:3,col)', result1(11:12,col)', [0, 0]...
result1(8:10,col), result1(4:7,col)', [20, 20]});

s2 = SuperQuadrics({result2(1:3,col)', result2(11:12,col)', [0, 0]...
   result2(8:10,col), result2(4:7,col)',[20,20]});

% exactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, 1e+03)
visualize_bounding_ellip(s1, s2);

% [prob, time] = enlargedBoundingVolume(s1, s2, Sigma1, Sigma2, 0.99, true)

visualizePositionErrors(s1, s2, Sigma1, Sigma2)

% 
mx = s2.tc-s1.tc;
[p_center, t_center, a_center, x_center] = linearChanceConstraintBound(s1, s2, mx, Sigma2, 'center-point', true)
% % [p_closed, t_closed, a_closed, x_closed] = linearChanceConstraintBound(s1, s2, mx, Sigma, 'closed-point', true)
% 
% visualize_bounding_ellip(s1, s2);
% [p_tangent, t_tangent, a_tangent, x_tangent] = linearChanceConstraintBound(s1, s2, mx, Sigma, 'tangent-point', true)
% 
% [p_tangent, t_tangent, a_tangent, x_tangent, r_tangent_ellip] = linearChanceConstraintBound(s1, s2, mx, Sigma, 'tangent-point-ellip', true)
% 
% [R,~,~] = svd(Sigma);
% s_x = SuperQuadrics({[r_tangent_ellip*Sigma(1,1), r_tangent_ellip*Sigma(2,2), r_tangent_ellip*Sigma(3,3)], [1,1], [0,0],...
%     s2.tc, rotm2quat(R), [20,20]});
% 
% s_x.PlotShape('r', 0.5)
% 
% [p_SQ, t_SQ, pdf_SQ, x_max] = maxPDFSQ(s1, s2, mx, Sigma)

% [p_sphere, t_sphere, pdf_sphere, x_max] = maxPDFSphere(s1, s2, mx, Sigma)

% scatter3(x_max(1), x_max(2), x_max(3), 'MarkerFaceColor', 'g');

% visualize_position_error(s1, s2, Sigma)

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