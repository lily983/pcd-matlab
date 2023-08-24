function prob = max_prob_gaussian_mixture(s1, s2, mx, Sigmax, method1, method2)
% Check dimension 
dimension = size(Sigmax, 2);

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

if strcmp(s1objectType, 'superquadric')==1 && strcmp(s2objectType, 'superquadric')==1
    prob = NaN;
    error('Input objects are not sphere or ellip, unable to use PCD-GMM');
end

%choose which bounding ellipsoid methods to use
switch method1
    case 'fast-computation'
        [mf, Sigmaf] = get_bounding_ellip(s1, s2);
    case 'accurate-result'
        [mf, Sigmaf] = get_bounding_ellip_fixed_point(s1, s2);
end

%choose which gmm to use(2 SG or 3SG or 4SG or 5SG)
[a, b] = get_gmm_param(Sigmaf, method2);

% compute probability of collision 
prob=0;

for i=1:size(a,2)
    prob = prob + a(i)*Gaussian_integration(mf, Sigmaf/b(i), mx, Sigmax) / b(i)^(dimension/2);
end

end
