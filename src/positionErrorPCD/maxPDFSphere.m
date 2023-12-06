function   [prob, time, pdf_max_x, x_max] = maxPDFSphere(s1, s2, mx, Sigma)
% maxPDFSphere: implementation of paper "Fast and Bounded Probabilistic
% Collision Detection for High-DOF Trajectory Planning in Dynamic Environments"
%
%Warning: this function can only be used for two sphere objects (because
%the searching region is the sum of two spheres
%
%Inputs
%   mx: Mean of relative position error x = x2-x1 (m = s2.tc-s1.tc)
%   Sigma: Covariance of the probability
%   s1, s2: Two convex bodies
%Outputs
%   x_max: Center position of s2 that has the maximum probability of
%   colliding with s1
%   pdf_max_x: PDF value at x_max
%   prob: The probability 
%   time: computation time

%Check if s1 and s2 are sphere 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);
 if strcmp(s1objectType, 'sphere')==0 && strcmp(s2objectType, 'sphere')==0
    pdf_max_x = NaN;
    error('Input objects are not sphere, unable to use Maxpdf');
 end

 tic;
 prob  = 0;

if collision_cfc(s1, s2)
    x_max = s2.tc;
    pdf_max_x = mvnpdf(x_max, mx, Sigma);
    return
end

opt =  optimoptions("fsolve","OptimalityTolerance",1e-15);
fun = @(lamda)norm((inv(Sigma) + lamda*eye(3))\(Sigma\mx + lamda*s1.tc) - s1.tc) - s1.a(1) - s2.a(1);
lamda0 = 1;
lamda = fsolve(fun, lamda0, opt);

x_max =double((inv(Sigma) + lamda*eye(3))\(Sigma\mx + lamda*s1.tc));
pdf_max_x = mvnpdf(x_max, mx, Sigma);

sphereVolume =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;

prob = sphereVolume * pdf_max_x;
time = toc;
end