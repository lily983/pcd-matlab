function x_new = method7(obj, R)
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
    results(i) = a - 1/norm(A\R(:,:,i)*ei);
end

[d, ~] = max(results);

x_new = (A^2 *n) ./ vecnorm(A * n, 2, 1) + d .* n;

end