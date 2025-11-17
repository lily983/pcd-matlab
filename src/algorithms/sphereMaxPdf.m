function   [prob, time, pdf_max_x, x_max] = sphereMaxPdf(e1, e2, xx, Sigma)
% sphereMaxPdf: implementation of paper "Fast and Bounded Probabilistic
% Collision Detection for High-DOF Trajectory Planning in Dynamic Environments"
%
%Warning: this function can only be used for two sphere objects (because
%the searching region is the sum of two spheres
%
%
%Inputs
%   e1, e1: Two soherical objects
%   xx: Mean of the position error x2. Noticed that this paper only support
%   the second object subjected to position error
%   Sigmax: Covariance of the relative position error
%Outputs
%   prob: The probability approximation
%   time: computation time


%Check if e1 and e1 are sphere 
e1objectType = getObjectType(e1);
e2objectType = getObjectType(e2);
 if strcmp(e1objectType, 'sphere')==0 && strcmp(e2objectType, 'sphere')==0
    pdf_max_x = NaN;
    error('Input objects are not sphere, unable to use Maxpdf');
 end

 tic;
 prob  = 0;

if collision_cfc(e1, e2)
    x_max = e2.tc;
    pdf_max_x = mvnpdf(x_max, xx, Sigma);
    return
end

opt =  optimoptions("fsolve","OptimalityTolerance",1e-15);

% The objective function is to find the surface point on the Minkowski sum S1+S2 centered at s1.tc which has the
% maximum pdf value at N(s2.tc, Sigmax)
fun = @(lamda)norm((inv(Sigma) + lamda*eye(3))\(Sigma\xx + lamda*e1.tc) - e1.tc) - e1.a(1) - e2.a(1);
lamda0 = 1;
lamda = fsolve(fun, lamda0, opt);

x_max =double((inv(Sigma) + lamda*eye(3))\(Sigma\xx + lamda*e1.tc));
pdf_max_x = mvnpdf(x_max, xx, Sigma);

sphereVolume =  4*pi/3*(e1.a(1)+e2.a(1))^3 ;

prob = sphereVolume * pdf_max_x;
time = toc;
end