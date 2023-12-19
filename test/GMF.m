function y = GMF(x, mu, Sigma, a, b, n)
% Given parameters of GMF function (ma, a, b) and dimension of space (n),
% and inputs x, return GMF results y
% Inputs:
% x: inputs, Nx1 vector, N is the size of x
% a, b: GMF parameters
% n: dimension of space
% mu, Sigma: center and defining matrix of the ellipsoid. If 1D, both are
% single numbers.
% Outputs:
% y: outpus, 1xN vector
sizeParameters = max(size(a));

y = zeros(1, max(size(x)));

for i = 1:sizeParameters
    for j = 1:max(size(x))
        y(j) = y(j) + a(i) * b(i)^(n/2) * exp( (-b(i)/2) * (x(j)-mu)' / Sigma * (x(j)-mu));
    end
end

y = y ./ ( (2*pi)^(n/2) * det(Sigma)^0.5);

end