% Define parameters for the convolved Gaussian (from f and h):
mu0 = [0, 0];         % Example mean
Sigma0 = [1, 0.3; 
          0.3, 2];   % Example covariance

% Define parameters for the ellipse:
c = [0, 0];   % Center of the ellipse
a = 2;        % Semi-axis length in x-direction
b = 3;        % Semi-axis length in y-direction

% Define the Gaussian PDF as an anonymous function:
gauss_pdf = @(x,y) (1/(2*pi*sqrt(det(Sigma0)))) * ...
    exp(-0.5 * ([x(:) y(:)] - mu0) * (Sigma0\([x(:) y(:)] - mu0))');

% For vectorized evaluation, we re-define it so that it works elementwise:
gauss_pdf_vec = @(x,y) (1/(2*pi*sqrt(det(Sigma0)))) * ...
    exp(-0.5*((x-mu0(1)).^2*Sigma0(2,2) - 2*(x-mu0(1)).*(y-mu0(2))*Sigma0(1,2) + (y-mu0(2)).^2*Sigma0(1,1)) ./ det(Sigma0));

% Define the indicator function for the ellipse:
indicator = @(x,y) (( (x - c(1)).^2 / a^2 + (y - c(2)).^2 / b^2 ) <= 1);

% Define the integration region.
% Since the ellipse is bounded, we can integrate over its bounding box:
x_min = c(1) - a;
x_max = c(1) + a;
y_min = c(2) - b;
y_max = c(2) + b;

% Define the integrand for normalization constant Z:
integrand_Z = @(x,y) gauss_pdf_vec(x,y) .* indicator(x,y);

% Compute Z (the probability mass inside the ellipse)
Z = integral2(integrand_Z, x_min, x_max, y_min, y_max, 'Method', 'iterated');

% Define integrands for the first moments:
integrand_x = @(x,y) x .* gauss_pdf_vec(x,y) .* indicator(x,y);
integrand_y = @(x,y) y .* gauss_pdf_vec(x,y) .* indicator(x,y);

% Compute the mean (first moment)
mu_x = integral2(integrand_x, x_min, x_max, y_min, y_max, 'Method', 'iterated') / Z;
mu_y = integral2(integrand_y, x_min, x_max, y_min, y_max, 'Method', 'iterated') / Z;
mu_trunc = [mu_x, mu_y];

% Define integrands for the second moments (for covariance)
integrand_xx = @(x,y) (x - mu_x).^2 .* gauss_pdf_vec(x,y) .* indicator(x,y);
integrand_yy = @(x,y) (y - mu_y).^2 .* gauss_pdf_vec(x,y) .* indicator(x,y);
integrand_xy = @(x,y) (x - mu_x).*(y - mu_y) .* gauss_pdf_vec(x,y) .* indicator(x,y);

cov_xx = integral2(integrand_xx, x_min, x_max, y_min, y_max, 'Method', 'iterated') / Z;
cov_yy = integral2(integrand_yy, x_min, x_max, y_min, y_max, 'Method', 'iterated') / Z;
cov_xy = integral2(integrand_xy, x_min, x_max, y_min, y_max, 'Method', 'iterated') / Z;
Sigma_trunc = [cov_xx, cov_xy; cov_xy, cov_yy];

% Display the results:
disp('Approximated (truncated) mean:');
disp(mu_trunc);
disp('Approximated (truncated) covariance:');
disp(Sigma_trunc);
