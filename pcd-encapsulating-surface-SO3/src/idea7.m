function x_new = idea7(obj, R, scale)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

% x_original  = obj.GetPoints();

x_new = zeros(3, obj.N(1) * obj.N(2));

% R_new = zeros(3,3);

n = obj.GetNormals();

for i = 1:size(R,3)
%     R_new = R_new + R(:,:,i);
x_new = x_new + R(:,:,i) * obj.GetPointsFromNormal(R(:,:,i) *n);
end

x_new = x_new  ./ scale;
    
end