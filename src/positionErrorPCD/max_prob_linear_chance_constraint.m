function prob = max_prob_linear_chance_constraint(s1, s2, Sigmax)
% Get bounding ellipsoid for the Minkowski sum of s1 and s2
[~, Sigmaf] = get_bounding_ellip(s1, s2);

% In this setting, we let mx = m2-m1 and mf = 0
mx = s2.tc-s1.tc;

% Transform the defining matrix of the bounding ellipsoid to identity
% (spheres at the origin)
mx_transform = Sigmaf^0.5 \ mx;

prob = 1/2 + 1/2*erf( (1-mx_transform'*

end