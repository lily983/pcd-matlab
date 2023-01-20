function prob = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2)
%max_contact_probability_polyhedron_pure_translation Compute the upper bound of contact
%probability for two polyhedrons objects
%
%Inputs:
%Sigma: The convariance matrix of the probability distribution for the
%position error of s2
%
%Outputs:
%prob: The upper bound of contact probability

prob = 0;

m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

% Scale and shift all points by 'Sigma \ (x-s2.tc)' because in this case
% the position error is for s2 not s1

mink_points = (sqrtm(Sigma)\x_mink)';
shift = sqrtm(Sigma) \ s2.tc;

mink_points(:,1) = mink_points(:,1) - shift(1,1);
mink_points(:,2) = mink_points(:,2) - shift(2,1);
mink_points(:,3) = mink_points(:,3) - shift(3,1);

X_ = reshape(mink_points(:,1), [20,20]); 
Y_= reshape(mink_points(:,2),  [20,20]); 
Z_ = reshape(mink_points(:,3),  [20,20]); 

patch_mink = surf2patch(X_,Y_,Z_, 'triangles');
%patch(patch_mink, 'FaceAlpha', 0.3)

s1_points = (sqrtm(Sigma) \ s1.GetPoints())';
s1_points(:,1) = s1_points(:,1) - shift(1,1);
s1_points(:,2) = s1_points(:,2) - shift(2,1);
s1_points(:,3) = s1_points(:,3) - shift(3,1);

s2_points = (sqrtm(Sigma) \ s2.GetPoints())';
s2_points(:,1) = s2_points(:,1) - shift(1,1);
s2_points(:,2) = s2_points(:,2) - shift(2,1);
s2_points(:,3) = s2_points(:,3) - shift(3,1);

patch_s1 = surf2patch(reshape(s1_points(:,1), [20,20]), reshape(s1_points(:,2), [20,20]), reshape(s1_points(:,3), [20,20]), 'triangles');
patch_s2 = surf2patch(reshape(s2_points(:,1), [20,20]), reshape(s2_points(:,2), [20,20]), reshape(s2_points(:,3), [20,20]), 'triangles');

[dist,pts,G,H] = GJK_dist(patch_s1, patch_s2);

n_d = (G - H)/norm(G-H);

prob = 0;

for i = 1:size(patch_mink.faces,1)
    v1 = mink_points(patch_mink.faces(i,1),:)';
    v2 = mink_points(patch_mink.faces(i,2),:)';
    v3 = mink_points(patch_mink.faces(i,3),:)';
    
    m = -cross(v1-v2, v1-v3);
    n = m/norm(m);
    area = 1/2*norm(m);
    
    %Compute max F(vi,n_d)*n
    F_array = zeros(1,3);
    F_array(1) = dot(F(v1,n_d),n);
    F_array(2) = dot(F(v2,n_d),n);
    F_array(3) = dot(F(v3,n_d),n);
    if isnan(max(F_array)*area)
        warning('find nan')
        v1
        v2
        v3
        continue
    end
    prob = prob + max(F_array)*area;
    
end

end

function result = F(x,n_d)

result = 1/(2*pi)*(1+erf(dot(x,n_d)/sqrt(2)))*n_d;

end