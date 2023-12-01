function [prob, t] = max_prob_ellip_UB_exact(s1, s2, Sigma)
% This function reproduce methods using CDF of quadratic form of Gaussian
% distributed variable to get PCD value. The code is writen based on the
% equations on paper.
% Difference between this function with quadraticExact() is that the
% condition to break the iteration and the place of (-1)^k
tic;
prob = 0;
    
mu = s2.tc-s1.tc;

% In Thomas paper, only spheres are considered P = P(x' * x < (r1+r2)^2)
% In Pan Jia paper, ellipsoids are considered P = P(x' * A * x < 1)
% This code does a global space transformation so that A = eye(3), r1+r2=1
% In this case, we will do y = inv(A^1/2) x to the
% space
[~, A] = get_bounding_ellip(s1, s2);

mu_transformed = A^0.5 \ mu;
Sigma_transformed = A^0.5 \ Sigma / A^0.5;
A_transformed = A^0.5 \ A / A^0.5; % eye(3)

[V,D] = eig(Sigma_transformed^1/2 * A_transformed * Sigma_transformed^1/2);

lambda = diag(D);

b = V' / Sigma_transformed^0.5 * mu_transformed;

prob = fixed_point(b, diag(lambda));
if imag(prob)>1e-02
    prob = NaN;
else
    prob = real(prob);
end

t = toc;
end


function prob = fixed_point(b, lamda)
prob=0;

% matlab array starts from position 1 instead of 0. So k starts from 1
% while in paper it starts from position 0
k=1;
c0 = exp(-0.5 * (b(1)^2 + b(2)^2 + b(3)^2)) * ((2*lamda(1)) * (2*lamda(2)) * (2*lamda(3)))^(-0.5);
c_list(k) = c0;

delta = (-1)^(k-1)*c_list(k)/gamma(3/2+k);
prob = prob+delta;

k=2;
c_next = fixed_point_iteration(c_list, b, k, lamda);

% keep getting more series until the increment is less than threshold
while k<30 && abs(delta)<1e-02
    c_list(k) = c_next; 
    delta = (-1)^(k-1)*c_list(k)/gamma(3/2+k);
    prob = prob+delta;
    
    k=k+1;
    c_next = fixed_point_iteration(c_list, b, k, lamda);
end

end

function c_k = fixed_point_iteration(c, b, k, lamda)
c_k=0;

for i=0:1:(k-2)
    c_k = c_k + get_dk(b, k-i, lamda) * c(i+1);
end

c_k = (1.0/(k-1)) * c_k;
end

function d_k = get_dk(b, k, lamda)
d_k = 0;

for i=1:3
    d_k = d_k + (1 - (k-1)*b(i)^2)*(2*lamda(i))^(-(k-1));
end

d_k = d_k * 0.5;
end