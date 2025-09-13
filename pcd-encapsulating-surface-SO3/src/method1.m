function x_new = method1(obj, R)
%idea: 1, 2, 3, 4
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

x_new = zeros(3, obj.N(1) * obj.N(2));

m = obj.GetGradients();

u = obj.GetHypersphereFromGradient(m);

A = diag(obj.a);

sum_RA = zeros(3);

for i = 1:size(R,3)
    sum_RA = sum_RA + R(:,:,i) * A;
end

x_new = sum_RA * u;

x_new = x_new ./ size(R,3);
    
end