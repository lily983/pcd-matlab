function prob = max_prob_gaussian_mixture(s1, s2, mx, Sigmax, method1)
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

[a1, b1] = get_gmm_param(Sigmaf, '1SG');

% compute probability of collision 
prob1=0;

for i=1:size(a1,2)
    prob1 = prob1 + a1(i)*Gaussian_integration(mf, Sigmaf/b1(i), mx, Sigmax) / b1(i)^(dimension/2);
end

[a2, b2] = get_gmm_param(Sigmaf, '5SG');

% compute probability of collision 
prob2=0;

for i=1:size(a2,2)
    prob2 = prob2 + a2(i)*Gaussian_integration(mf, Sigmaf/b2(i), mx, Sigmax) / b2(i)^(dimension/2);
end

if prob1>prob2
    prob=prob2;
else
    prob=prob1;
end

end
