function x_new = method5(obj, R, k)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m
x_new = zeros(3, obj.N(1) * obj.N(2));

sum_R = zeros(3);

for i = 1:size(R,3)
    sum_R = sum_R + R(:,:,i);
end
sq_eye3 = SuperQuadrics({obj.a, obj.eps, [0,0], zeros(3,1), rotm2quat(eye(3)), obj.N});

x_new = sum_R * sq_eye3.GetPoints();

x_new = k * x_new ./ size(R,3);
end