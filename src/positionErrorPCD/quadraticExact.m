function [prob, t] = quadraticExact(s1, s2, mx, Sigma)
% This function reproduce methods using CDF of quadratic form of Gaussian
% distributed variable to get PCD value. The code is modified based on
% source code provided by Anatony Thomas https://bitbucket.org/1729antony/ comparison_cp_methods/src/master/
    % Inputs:
        % s1, s2: sphere or ellipsoid
        % mx: center of position error distribution
        % Sigma: covariance matrix of position error distribution
    % Outputs:
        % prob: PCD value
        % t: computation time
tic;
prob = 0;

% In Thomas paper, only spheres are considered P = P(x' * x < (r1+r2)^2)
% In Pan Jia paper, ellipsoids are considered P = P(x' * A * x < 1)
% This code does a global space transformation so that A = eye(3), r1+r2=1
% In this case, we will do y = inv(A^1/2) x to the
% space
radius = 1;
[~, A] = get_bounding_ellip(s1, s2);

mu_transformed = A^0.5 \ mx;
Sigma_transformed = A^0.5 \ Sigma / A^0.5;
A_transformed = A^0.5 \ A / A^0.5; % eye(3)

[V,D] = eig(Sigma_transformed^1/2 * A_transformed * Sigma_transformed^1/2);

lambda = diag(D);

b = V' / Sigma_transformed^0.5 * mu_transformed;

j = 0;
tmp = 0;

while (true)

    prob = prob + c_param(b,j,lambda)* (radius^(1.5+j))/gamma(j+2.5);
    j = j +1;

    % to avoid the anomaly, which might happen when sigma very very negligible 
    if prob < 1 && prob >= 0
        tmp = prob;
    end
    if prob < 1e-50
        break;
    end
    if j > 15
        prob = tmp;
        break;
    end

end   

prob = tmp;
t = toc;
end


function value = d_param(b, k, lambda)
value = (-1)^k *0.5 * ((1-k*b(1)^2)*(2*lambda(1))^(-k) +  (1-k*b(2)^2)*(2*lambda(2))^(-k) + (1-k*b(3)^2)*(2*lambda(3))^(-k));
end

function value = c_param(b,k,lambda)
value = 0;
if (k == 0)
    value = exp(-0.5*(b(1)^2 + b(2)^2 + b(3)^2))*(2*lambda(1))^(-0.5)*(2*lambda(2))^(-0.5)*(2*lambda(3))^(-0.5);
else
    for i = 0:k-1
        value = value + d_param(b, k-i, lambda)*c_param(b,i, lambda);
    end
    value = value/k;
end
end
