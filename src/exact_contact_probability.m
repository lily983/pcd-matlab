function prob = exact_contact_probability(mu, Sigma, s1, s2, N)
%   exact_contact_probability: Numerical approximation of
%   the exact probability of collision
%   
%   Inputs:
%   mu: Mean pose (center of s2)
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%   N: Sampling points number
prob = 0;

mu_log_skew = logm(mu);
mu_log = [vex(mu_log_skew(1:3,1:3)); mu_log_skew(1:3,4)];

samples = mvnrnd(mu_log', Sigma, N);

s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

for i=1:N
    samples_skew = [skew(samples(i,1:3)), samples(i,4:6)'; zeros(1,4)];
    g_rand = expm(samples_skew);
    
    s3.tc = g_rand(1:3,4);
   % for two sphere objects, use simplier method to test if they collide
    if s1.eps(1)==1 && s1.eps(2)==1 && s3.eps(1) ==1 && s3.eps(2)==1 && isequal(s1.a./s1.a(1), ones(1,3))==true && isequal(s3.a./s3.a(1), ones(1,3))==true
       if norm(s1.tc-s3.tc)<=(s1.a(1)+s3.a(1))
           prob = prob+1;      
       end
    else
       if collision_cfc(s1,s3,'constrained')
           prob = prob+1;
       end
   end
end
prob=prob/N;

end