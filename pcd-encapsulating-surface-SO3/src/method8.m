function x_new = method8(obj, R)
%%%%Ellipsoid only%%%%
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m
n = obj.GetNormals();

A = diag(obj.a);

[~, idx] = max(obj.a);

ei = zeros(3,1);
ei(idx) = 1;

results = zeros(1, size(R,3));

for i = 1:size(R,3)
    Ri = R(:,:,i);
    results(i) = (ei' * Ri' / A^2 * Ri * A * ei) / vecnorm(A^2 \ Ri * ei, 2, 1) ...
        - (ei' * Ri' / A^2 * Ri * ei) / (vecnorm(A^2 \ Ri * ei, 2, 1) * vecnorm(A^1 \ Ri * ei, 2, 1)) ;
end

[p, ~] = max(results);

x_new = (A^2 *n) ./ vecnorm(A * n, 2, 1) + p .* n;

end