function prob = max_prob_gaussian_mixture(s1, s2, mx, Sigmax, aRatio, bRatio, mb)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

[a1,a2,b1,b2] = get_gmm_param(Sigmaf, aRatio, bRatio, mb);

dimension = size(Sigmax,1);

prob = a1 * Gaussian_integration(mf, Sigmaf/b1, mx, Sigmax) / b1^(dimension/2) ...
    + a2 * Gaussian_integration(mf, Sigmaf/b2, mx, Sigmax) / b2^(dimension/2) ;
end