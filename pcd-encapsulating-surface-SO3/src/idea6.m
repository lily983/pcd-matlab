function x_new = idea6(obj, R, mu)
% weighted sum; weight depends on the deviation from the mean rotation
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

x_original  = obj.GetPoints();

x_new = zeros(3, obj.N(1) * obj.N(2));

R_new = zeros(3,3);

for i = 1:size(R,3)
    d(i) = norm(skew2vec(logm_SO(R(:,:,i)' * mu)));
end
D = sum(d);

for i = 1:size(R,3)
    R_new = R_new + d(i) * R(:,:,i);
end

R_new = R_new ./ D;

x_new = R_new * x_original;
    
end