function   [pdf_max_x, x_max] = max_probability_pure_translation(mu, Sigma, s1, s2)
% max_probability_pure_translation: Center position of s2 is a Guassian random variable, 
% the function calculate the maximum probability that it will collide with s1 
%
%Inputs
%   mu: Mean pose (center of s2)
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%Outputs
%   x_max: Center position of s2 that has the maximum probability of
%   colliding with s1
%   pdf_max_x: The probability of x_max

if collision_cfc(s1,s2,'constrained')
    x_max = s2.tc;
    pdf_max_x  = mvnpdf(x_max, mu, Sigma);
    return
end

syms lamda
eqn = norm((inv(Sigma) + lamda*eye(3))\(Sigma\mu + lamda*s1.tc) - s1.tc) == s1.a(1) + s2.a(1);
sol = solve(eqn,lamda,'Real',true);
sol = double(sol(2,1));

x_max = (inv(Sigma) + sol*eye(3))\(Sigma\mu + sol*s1.tc);

pdf_max_x = mvnpdf(x_max, mu, Sigma);
pdf_max_x = double(pdf_max_x);
end