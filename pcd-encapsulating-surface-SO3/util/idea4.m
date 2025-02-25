function x_new = idea4(obj, R)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

x_original  = obj.GetPoints();

x_new = zeros(3, obj.N(1) * obj.N(2));

R_new = zeros(3,3);

for i = 1:size(R,3)
    R_new = R_new + R(:,:,i);
end

R_new = R_new ./ (size(R,3)-2);

x_new = R_new * x_original;
    
end