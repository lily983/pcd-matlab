function x_new = method10(obj, R)
%idea: 1, 2, 3, 4
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

x_new = zeros(3, obj.N(1) * obj.N(2));

n = obj.GetNormals();

for i = 1:size(R,3)
x_new = x_new + R(:,:,i) * obj.GetPointsFromNormal(R(:,:,i)' *n);
end

x_new = x_new  ./ scale;
end