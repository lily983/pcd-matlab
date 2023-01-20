function prob = exact_contact_probability_pure_translation(mu, Sigma, s1, s2, N)
%   exact_contact_probability_pure_translation: Numerical approximation of
%   the exact probability of collision
%   
%   Inputs:
%   mu: Center position of s2
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%   N: Sampling points number
prob = 0;

samples = mvnrnd(mu, Sigma, N);

s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

for i=1:N
    
    s3.tc  = samples(i,:)';
    
    % if s1 and s2 are sphere. If they are, using simplier method to
    % compute the collision status
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



