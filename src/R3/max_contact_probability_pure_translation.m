function   [pdf_max_x, x_max] = max_contact_probability_pure_translation(mu, Sigma, s1, s2)
% max_probability_pure_translation: Center position of s2 is a Guassian random variable, 
% the function calculate the maximum probability that it will collide with s1 
%
%Warning: this function can only be used for two sphere objects (because
%the searching region is the sum of two spheres
%
%Inputs
%   mu: Mean pose (center of s2)
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%Outputs
%   x_max: Center position of s2 that has the maximum probability of
%   colliding with s1
%   pdf_max_x: The probability of x_max

%Check if s1 and s2 are sphere 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);
 if strcmp(s1objectType, 'sphere')==0 && strcmp(s2objectType, 'sphere')==0
    pdf_max_x = NaN;
    error('Input objects are not sphere, unable to use PCD-maxpdf');
end

if collision_cfc(s1, s2)
    x_max = s2.tc;
    pdf_max_x = mvnpdf(x_max, mu, Sigma);
    return
end

opt =  optimoptions("fsolve","OptimalityTolerance",1e-15);
fun = @(lamda)norm((inv(Sigma) + lamda*eye(3))\(Sigma\mu + lamda*s1.tc) - s1.tc) - s1.a(1) - s2.a(1);
lamda0 = 1;
lamda = fsolve(fun, lamda0, opt);

x_max =double((inv(Sigma) + lamda*eye(3))\(Sigma\mu + lamda*s1.tc));
pdf_max_x = mvnpdf(x_max, mu, Sigma);

end