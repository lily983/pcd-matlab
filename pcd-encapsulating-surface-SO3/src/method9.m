function x_new = method9(obj, R, k)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

x_new = zeros(3, obj.N(1) * obj.N(2));

n = obj.GetNormals();

for i = 1:size(R,3)
x_new = x_new + R(:,:,i) * obj.GetPointsFromNormal(R(:,:,i) *n);
end

x_new = k * x_new  ./ size(R,3);
end