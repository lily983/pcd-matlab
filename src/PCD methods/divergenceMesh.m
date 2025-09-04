function [prob, time] = divergenceMesh(s1, s2, mx, Sigma)
% divergenceMesh: implementation of paper "Efficient Probabilistic Collision Detection for Non-Convex Shapes"  
%
%Inputs
%   mx: Mean of relative position error x = x2-x1 
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%Outputs
%   prob: The probability 
%   time: computation time
% Check surface sampling points
if s1.N ~= s2.N
    prob = NaN;
    error("S1 and S2 sampling points (N) are not the same");
end

% tic;
% prob = 0;

% Surface points
SN = s1.N;

m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
% S1 and S2 are subjected position errors x1~N(m1, Sigma1) and x2~N(m2,
% Sigma2). Here S1 and S2 are assumed to center at origin
x_mink = mink.GetMinkSumFromGradient(m1);

% Scale and shift all points by 'sqrtm(Sigma) \ (x-mx)' so that position
%error center at origin and Sigma is eye(3), a unit ball
mink_points = (sqrtm(Sigma)\x_mink)';
shift = sqrtm(Sigma) \ mx;

mink_points(:,1) = mink_points(:,1) - shift(1,1);
mink_points(:,2) = mink_points(:,2) - shift(2,1);
mink_points(:,3) = mink_points(:,3) - shift(3,1);

X_ = reshape(mink_points(:,1), SN); 
Y_= reshape(mink_points(:,2),  SN); 
Z_ = reshape(mink_points(:,3),  SN); 

patch_mink = surf2patch(X_,Y_,Z_, 'triangles');
%patch(patch_mink, 'FaceAlpha', 0.3)

tic;
prob = 0;


% In this paper, direction vector n_d is the connection between closed
% points from Sigma^0.5*s2 - Sigma^0.5*s1
s1_points = (sqrtm(Sigma) \ s1.GetPoints())';
s2_points = (sqrtm(Sigma) \ s2.GetPoints())';

patch_s1 = surf2patch(reshape(s1_points(:,1), SN), reshape(s1_points(:,2), SN), reshape(s1_points(:,3), SN), 'triangles');
patch_s2 = surf2patch(reshape(s2_points(:,1), SN), reshape(s2_points(:,2), SN), reshape(s2_points(:,3), SN), 'triangles');

[~,~,G,H] = GJK_dist(patch_s1, patch_s2);

n_d = (G - H)/norm(G-H);

% prob = 0;

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
        
        continue
    end
    prob = prob + max(F_array)*area;
    
end
time = toc;
end

function result = F(x,n_d)

result = 1/(2*pi)*(1+erf(dot(x,n_d)/sqrt(2)))*n_d;

end