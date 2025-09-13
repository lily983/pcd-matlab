function x_new = method5(obj, R)
%idea: 1, 2, 3, 4
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m
x_new = zeros(3, obj.N(1) * obj.N(2));

n = obj.GetNormals();

A = diag(obj.a);

sum_R = zeros(3);

for i = 1:size(R,3)
    sum_R = sum_R + R(:,:,i);
end

x_new = x_new ./ size(R,3);
end