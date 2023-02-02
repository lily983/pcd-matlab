function prob = exact_contact_probability(mu, Sigma, s1, s2, N)
%   exact_contact_probability: Numerical approximation of
%   the exact probability of collision
%   
%   Inputs:
%   mu: Mean pose (center of s2) 4by4 matrix
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%   N: Sampling points number
prob = 0;

mu_vee = get_vee_vector(mu);

samples = mvnrnd(mu_vee', Sigma, N);

s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

for i=1:N
    g_rand = get_SE3_matrix(samples(i, 1:6)');
    
    s3.tc = g_rand(1:3,4);
   % for two sphere objects, use simplier method to test if they collide
    if s1.eps(1)==1 && s1.eps(2)==1 && s3.eps(1) ==1 && s3.eps(2)==1 && isequal(s1.a./s1.a(1), ones(1,3))==true && isequal(s3.a./s3.a(1), ones(1,3))==true
       if norm(s1.tc-s3.tc)<=(s1.a(1)+s3.a(1))
           prob = prob+1;      
       end
    else
       if collision_cfc(s1,s3)
           prob = prob+1;
       end
   end
end
prob=prob/N;

end