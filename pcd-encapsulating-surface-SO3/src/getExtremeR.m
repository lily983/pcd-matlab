function [test, dist] = getExtremeR(R, mean,Sigma)
dist = zeros(1, size(R,3));

for i = 1:size(R,3)
    dist(i) =  norm(skew2vec(logm_SO(R(:,:,i)' * mean)));
end

chi2_val = chi2inv(0.95, 3); % Chi-squared value for 99% confidence interval in 3D

% Perform eigenvalue decomposition of the covariance matrix
[~, D] = eig(Sigma);

% The square roots of the eigenvalues give the axes lengths of the ellipsoid
radii = sqrt(diag(D)) *  sqrt(chi2_val);  % Radii of the ellipsoid

radii_dist = norm(skew2vec(logm_SO(expm_SO(radii) * mean)));

Q3 = prctile(dist, 75);  % Third quartile (75th percentile)

indices = find(dist > Q3);

test = zeros(3,3,size(indices,2));

for i = 1:size(indices,2)
    test(:,:,i) = R(:,:,indices(i));
end
    

end