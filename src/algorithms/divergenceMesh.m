function [prob, time] = divergenceMesh(s1, s2, xx, Sigmax)
% divergenceMesh: implementation of paper "Efficient Probabilistic Collision Detection for Non-Convex Shapes"  
%
%Inputs
%   s1, s2: Two convex mesh objects (converts superquadrics to mesh in
%   code)
%   xx: Mean of the relative position error x = x2-x1 
%   Sigmax: Covariance of the relative position error
%Outputs
%   prob: The probability approximation
%   time: computation time

% The paper supports mesh objects, here we first convert superquadrics to
% mesh 
s1_points = (sqrtm(Sigmax) \ s1.GetPoints())';
s2_points = (sqrtm(Sigmax) \ s2.GetPoints())';

% Get the mesh of superquadrics S1 and S2
patch_s1 = surf2patch(reshape(s1_points(:,1), s1.N), reshape(s1_points(:,2), s1.N), reshape(s1_points(:,3), s1.N), 'triangles');
patch_s2 = surf2patch(reshape(s2_points(:,1), s2.N), reshape(s2_points(:,2), s2.N), reshape(s2_points(:,3), s2.N), 'triangles');

% Because surface points of s1 and s2 may not be the same, here we use all pairwise sums (s1.N*s2.N x 3)
[I,J] = ndgrid(1:size(s1_points,1), 1:size(s2_points,1));
mink_points=s1_points(I(:),:) - s2_points(J(:),:);

% Surface of Minkowski sum is the convex hull in 3D
K = convhulln(mink_points);

patch_mink = struct('Faces', K, 'Vertices', mink_points);

% Shifts the Minkowski sum to the origin
patch_mink.Vertices = patch_mink.Vertices - (sqrtm(Sigmax) \ xx)'; 

% Start record algorithm running time
tic;
prob = 0;

% In this paper, direction vector n_d is the connection between closed
% points from Sigmax^0.5*s2 - Sigmax^0.5*s1. 
[~,~,G,H] = GJK_dist(patch_s1, patch_s2);

n_d = (G - H)/norm(G-H);

for i = 1:size(patch_mink.Faces,1)
    v1 = mink_points(patch_mink.Faces(i,1),:)';
    v2 = mink_points(patch_mink.Faces(i,2),:)';
    v3 = mink_points(patch_mink.Faces(i,3),:)';
    
    m = -cross(v1-v2, v1-v3);
    n = m/norm(m);
    area = 1/2*norm(m);
    
    %Compute max F(vi,n_d)*n
    F_array = zeros(1,3);
    F_array(1) = dot(F(v1,n_d),n);
    F_array(2) = dot(F(v2,n_d),n);
    F_array(3) = dot(F(v3,n_d),n);
    if isnan(max(F_array)*area)
        continue
    end
    prob = prob + max(F_array)*area;
    
end
time = toc;
end

function result = F(x,n_d)

result = 1/(2*pi)*(1+erf(dot(x,n_d)/sqrt(2)))*n_d;

end