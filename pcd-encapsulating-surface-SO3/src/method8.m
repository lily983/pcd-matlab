function x_new = method8(obj, R)
%idea: 1, 2, 3, 4
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m
n = obj.GetNormals();

[~, idx] = max(obj.a);

ei = zeros(3,1);
ei(idx) = 1;

results = zeros(1, size(R,3));

for i = 1:size(R,3)
    R = R(:,:,i);
    results(i) = (ei' * R' / A^2 * R * A * ei) / vecnorm(A^2 \ R * ei, 2, 1) ...
        + (ei' * R' / A^2 * R * ei) / (vecnorm(A^2 \ R * ei, 2, 1) * vecnorm(A^1 \ R * ei, 2, 1)) ;
end

[d, ~] = max(results);

x_new = (A^2 *n) ./ vecnorm(A * n, 2, 1) + d .* n;

end