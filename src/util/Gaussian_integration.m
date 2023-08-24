function f = Gaussian_integration(mu1, Sigma1, mu2, Sigma2)
% Get the integration result of the product of two Gaussian functions over
% the entire space. It spoort 2D and 3D case.
    dimension = size(Sigma1, 2);
    f = (2*pi)^(-dimension/2)  * (det(Sigma1) * det(Sigma2) * det(inv(Sigma1) + inv(Sigma2)))^(-1/2)...
        * exp(1/2 * (mu1' / Sigma1 + mu2' / Sigma2) / (inv(Sigma1) + inv(Sigma2)) * (Sigma1' \ mu1 + Sigma2' \ mu2)...
        - 1/2 * (mu1' / Sigma1 * mu1 + mu2' / Sigma2 * mu2) );
end