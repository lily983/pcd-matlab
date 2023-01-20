function   [prob_x_max, x_max] = max_probability_pure_translation(mu, Sigma, s1, s2)
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
%   prob_x_max: The probability of x_max

if collision_cfc(s1,s2,'constrained')
    x_max = s2.tc;
    prob_x_max  = mvnpdf(x_max, mu, Sigma)
    return
end

syms lamda
eqn = norm((inv(Sigma) + lamda*eye(3))\(Sigma\mu + lamda*s1.tc) - s1.tc) == s1.a(1) + s2.a(1);
sol = solve(eqn,lamda,'Real',true);
solution = double(sol(2));


x_max = (inv(Sigma) + solution*eye(3))\(Sigma\mu + solution*s1.tc);

prob_x_max = mvnpdf(x_max, mu, Sigma);
prob_x_max = double(prob_x_max);
end