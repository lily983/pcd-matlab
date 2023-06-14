function [prob, weightedNorm] = max_prob_single_gaussian(s1, s2, mx, Sigmax)
% Check dimension 
dimension = size(s1.a, 2);

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

 if strcmp(s1objectType, 'superquadric')==1 && strcmp(s2objectType, 'superquadric')==1
    prob = NaN;
    weightedNorm = NaN;
    error('Input objects are not sphere or ellip, unable to use PCD-SG');
end

[mf, Sigmaf] = get_bounding_ellip(s1, s2);

k = (2*pi)^(dimension/2) * sqrt(det(Sigmaf) * exp(1));

% Gaussian intergration between bounding ellip and pdf
prob = k * Gaussian_integration(mf, Sigmaf, mx, Sigmax);
weightedNorm = mf'/Sigmaf*mf;
end