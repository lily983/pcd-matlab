function [prob, time] = gaussianMixturedFunction(s1, s2, mx, Sigmax, method)
% gaussianMixturedFunction: Closed-form solution for two ellipsoids
% subjected to position errors proposed by us.  
%Inputs
%   mx: Mean of relative position error x = x2-x1
%   Sigma: Covariance of the probability
%   s1, s2: SuperQuadrics
%Outputs
%   prob: Probability 
%   time: computation time

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

if strcmp(s1objectType, 'superquadric')==1 && strcmp(s2objectType, 'superquadric')==1
    error('Input objects are not sphere or ellip, unable to use PCD-GMM');
end

tic;
prob = 0;

if nargin == 4
    method = 'fast-computation';
end

%choose which bounding ellipsoid methods to use
switch method
    case 'fast-computation'
        [~, Sigmaf] = get_bounding_ellip(s1, s2);
    case 'accurate-result'
        [~, Sigmaf] = get_bounding_ellip_fixed_point(s1, s2);
end
mf = zeros(3,1);

[a1, b1] = get_gmm_param(Sigmaf, '2SG');

% compute probability of collision 
prob1=0;

for i=1:size(a1,2)
    prob1 = prob1 + a1(i)*Gaussian_integration(mf, Sigmaf/b1(i), mx, Sigmax) / b1(i)^(3/2);
end

[a2, b2] = get_gmm_param(Sigmaf, '5SG');

% compute probability of collision 
prob2=0;

for i=1:size(a2,2)
    prob2 = prob2 + a2(i)*Gaussian_integration(mf, Sigmaf/b2(i), mx, Sigmax) / b2(i)^(3/2);
end

if prob1>prob2
    prob=prob2;
else
    prob=prob1;
end

time = toc;
end
