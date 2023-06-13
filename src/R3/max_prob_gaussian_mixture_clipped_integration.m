function prob=max_prob_gaussian_mixture_clipped_integration(s1, s2, mx, Sigmax)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

[a1,a2,b1,b2] = get_gmm_param(Sigmaf, 10, 10, 3);

[mf, Sigmaf] = get_bounding_ellip(s1, s2);

% Check dimension 
dimension = size(Sigmax, 1);

kb1 = Gaussian_integration(mf, Sigmaf/b1, mx, Sigmax);
kb2 = Gaussian_integration(mf, Sigmaf/b2, mx, Sigmax);

prob = a1 * kb1* ellip_clipped_Gaussian_integration(1/sqrt(b1), dimension) / b1^(dimension/2) ...
    + a2 * kb2 * ellip_clipped_Gaussian_integration(1/sqrt(b2), dimension) / b2^(dimension/2) ;
end