function x_new = idea3(obj, R)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

x_new = zeros(3, obj.N(1) * obj.N(2));

n = obj.GetNormals();

A = diag(obj.a);

for i = 1:size(R,3)
    x_new =  x_new + (R(:,:,i) * A^2 * n)./ vecnorm(A * n) ;
end

x_new = x_new ./ size(R,3);
    
end