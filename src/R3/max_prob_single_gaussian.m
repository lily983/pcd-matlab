function [prob, weightedNorm] = max_prob_single_gaussian(s1, s2, mx, Sigmax)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

% Check dimension 
dimension = size(s1.a, 2);

k = (2*pi)^(dimension/2) * sqrt(det(Sigmaf) * exp(1));

% Gaussian intergration between bounding ellip and pdf
prob = k * Gaussian_integration(mf, Sigmaf, mx, Sigmax);
weightedNorm = mf'/Sigmaf*mf;
end