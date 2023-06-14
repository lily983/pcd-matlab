function prob = max_prob_gaussian_mixture(s1, s2, mx, Sigmax)
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

[a1,a2,b1,b2] = get_gmm_param(Sigmaf, 10, 3);

prob = a1 * Gaussian_integration(mf, Sigmaf/b1, mx, Sigmax) / b1^(dimension/2) ...
    + a2 * Gaussian_integration(mf, Sigmaf/b2, mx, Sigmax) / b2^(dimension/2) ;
end