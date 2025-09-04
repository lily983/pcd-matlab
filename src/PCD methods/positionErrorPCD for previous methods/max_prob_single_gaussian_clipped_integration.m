function prob=max_prob_single_gaussian_clipped_integration(s1, s2, mx, Sigmax)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

% Check dimension 
dimension = size(Sigmax, 1);

k = (2*pi)^(dimension/2) * sqrt(det(Sigmaf) * exp(1));

% Gaussian intergration between bounding ellip and pdf
kc = Gaussian_integration(mf, Sigmaf, mx, Sigmax);

invT = inv(Sigmaf)+inv(Sigmax);
% T=inv(invT);
% c = 1/(det(Sigmaf)*det(invT))^dimension;
% T = c.*inv(invT);

prob = k * kc * ellip_clipped_Gaussian_integration(sqrt(1/det(invT)), dimension);
end