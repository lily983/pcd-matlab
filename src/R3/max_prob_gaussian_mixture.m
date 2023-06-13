function prob = max_prob_gaussian_mixture(s1, s2, mx, Sigmax)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

[a1,a2,b1,b2] = get_gmm_param(Sigmaf, 10, 3);

dimension = size(Sigmax,1);

prob = a1 * Gaussian_integration(mf, Sigmaf/b1, mx, Sigmax) / b1^(dimension/2) ...
    + a2 * Gaussian_integration(mf, Sigmaf/b2, mx, Sigmax) / b2^(dimension/2) ;
end