close all;clear;clc;

% figure; hold on; axis equal; axis off;
figure; hold on; 

s1 = SuperEllipse([0.04, 0.2, 1, 0, 0.15, 0.35, 0.0, 50]);
s2 = SuperEllipse([0.05, 0.08, 1, 0, 0.18, 0.01, 0.3, 50]);

%%% Plot the effect of position errors
n = 20;

s1_color = hex2rgb('45AC59');
s1.PlotShape(s1_color, 1, 1);
Sigma1 = eye(2, 2) * 1 *1e-04;

x1_rand = mvnrnd(zeros(2,1), Sigma1, n);
s1_shift = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
    s1.tc(1), s1.tc(2), s1.ang, s1.N]);
for i = 1:size(x1_rand, 1)
    s1_shift.tc = s1.tc + x1_rand(i,:)';
    s1_shift.PlotShape(s1_color, 0.15, 0.0);
end

s2_color = hex2rgb('4FAAD1')
s2.PlotShape(s2_color, 1, 1);
Sigma2 =  eye(2, 2)  *1e-04;

x2_rand = mvnrnd(zeros(2,1), Sigma2, n);
s2_shift = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
    s2.tc(1), s2.tc(2), s2.ang, s2.N]);
for i = 1:size(x1_rand, 1)
    s2_shift.tc = s2.tc + x2_rand(i,:)';
    s2_shift.PlotShape(s2_color, 0.15, 0.0);
end

%% get upper bound ellipsoid
[mf, Sigmaf] = get_bounding_ellip(s1, s2);
[U, Lamda] = svd(Sigmaf);

s_up = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
            s2.tc(1), s2.tc(2), s2.ang, s2.N]);
s_up.a = sqrt(diag(Lamda));
s_up.ang = rotm2angle(U);
s_up.tc = mf;
s_up.N=100;

V_up = pi * det(Lamda)^0.5;

%% get the lower bound ellipse 1
[mf_lower_john, Sigmaf_lower_john] = get_lower_bounding_ellip(s1, s2);
[U_lower_john, Lamda_lower_john] = svd(Sigmaf_lower_john);

s_lower_john = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
            s2.tc(1), s2.tc(2), s2.ang, s2.N]);
s_lower_john.a = sqrt(diag(Lamda_lower_john));
s_lower_john.ang = rotm2angle(U_lower_john);
s_lower_john.tc = mf_lower_john;
s_lower_john.N=100;

V_lower_john  = pi * det(Lamda_lower_john)^0.5;

%%
figure; hold on

% Get minkowski sum of s1 and s2
s1Points = s1.GetPoints()' - s1.tc';
s2Points = -1*s2.GetPoints()' + s2.tc';
pgon1 = polyshape(s1Points(:,1), s1Points(:,2));
pgon2= polyshape(s2Points(:,1), s2Points(:,2));
mSum = minkowskiSum(pgon1, pgon2);

patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('45498C'), 'FaceAlpha', 0.0, 'EdgeAlpha', 1);
s_up.PlotShape(hex2rgb('70cf05'), 0.0, 0.8); % light green
s_lower_john.PlotShape(hex2rgb('ED5564'), 0.0, 0.8); % light red

%%% Put axes center at the origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.FontSize = 8;

%% compare the approximation and baseline
covApprox = mvnpdf([0 0], s1.tc' + s2.tc', Sigma1 + Sigma2 + Sigmaf)

p = (s2.tc-s1.tc);
a = p ./ norm(p);
Sigmax = Sigma1 + Sigma2;
prob_lcc = 0.5 + 0.5 * erf( (1 - norm(Sigmaf^0.5 \ p)) / sqrt(a' * (Sigmaf^0.5 \ Sigmax / Sigmaf^0.5) * a))

% [prob, time] = exactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, 1e+4)

% [prob, time] = fastExactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, 1e+5)

%% visualize the convolution and its Gaussian approximation
L = 2; % Total support for convolution from -L to L
a = L/3; % Support for each function from -a to a
N = 256; % Original grid size
M = 2^nextpow2(3*N-2); % Padded size, at least 3N-2

% Define original grid for each function from -a to a
x_original = linspace(-a, a, N);
[X_original, Y_original] = meshgrid(x_original, x_original);

% Vectorized computation of the Gaussian
X_shifted = X_original - s1.tc(1);
Y_shifted = Y_original - s1.tc(2);
det_Sigma = det(Sigma1); % Determinant of Sigma
inv_Sigma = inv(Sigma1); % Inverse of Sigma

% Compute the quadratic form (x - mu)' * Sigma^(-1) * (x - mu)
rho_1 = zeros(size(X_original));
for i = 1:N
    for j = 1:N
        x_vec = [X_shifted(i,j); Y_shifted(i,j)]; % Vector [x - mu_x, y - mu_y]
        rho_1(i,j) = exp(-0.5 * x_vec' * inv_Sigma * x_vec) / (2*pi * sqrt(det_Sigma));
    end
end

X_shifted = X_original + s2.tc(1);
Y_shifted = Y_original + s2.tc(2);
det_Sigma = det(Sigma2); % Determinant of Sigma
inv_Sigma = inv(Sigma2); % Inverse of Sigma

% Compute the quadratic form (x - mu)' * Sigma^(-1) * (x - mu)
rho_2 = zeros(size(X_original));
for i = 1:N
    for j = 1:N
        x_vec = [X_shifted(i,j); Y_shifted(i,j)]; % Vector [x - mu_x, y - mu_y]
        rho_2(i,j) = exp(-0.5 * x_vec' * inv_Sigma * x_vec) / (2*pi * sqrt(det_Sigma));
    end
end

inv_Sigma_f = inv(Sigmaf);
h = zeros(N, N);
for i = 1:N
    for j = 1:N
        x_vec = [X_original(i,j); Y_original(i,j)]; % [x, y] (centered at origin)
        quad_term = x_vec' * inv_Sigma_f * x_vec; % x' * Sigma_f^(-1) * x
        h(i,j) = double(quad_term <= 1); % 1 inside ellipse, 0 outside
    end
end

figure; hold on
surf(X_original, Y_original, h);
surf(X_original, Y_original, rho_1);
surf(X_original, Y_original, rho_2);


%%
p = (M - N)/2;
f_padded = padarray(rho_1, [p p], 'both');
g_padded = padarray(rho_2, [p p], 'both');
h_padded = padarray(h, [p p], 'both');

figure; hold on
% Plot the result
x_padded = linspace(-L, L, M);
[X_padded, Y_padded] = meshgrid(x_padded, x_padded);
surf(X_padded, Y_padded, f_padded);
surf(X_padded, Y_padded, g_padded);
surf(X_padded, Y_padded, h_padded);

%% Compute convolution
FFT_F = fft2(f_padded);
FFT_G = fft2(g_padded);
FFT_H = fft2(h_padded);
Y = FFT_F .* FFT_H .* FFT_G;
C = ifft2(Y);

% Plot the result
figure;
x_padded = linspace(-L, L, M);
[X_padded, Y_padded] = meshgrid(x_padded, x_padded);
surf(x_padded, x_padded, C);
xlabel('X'); ylabel('Y'); title('Convolution Result'); 
colorbar;
shading flat;

%% Plot the Gaussian approximation
mu_ap = -s2.tc+s1.tc;
Sigma_ap = Sigmaf + Sigma1 + Sigma2;

rho_ap = zeros(size(X_original));

X_shifted = X_original + mu_ap(1);
Y_shifted = Y_original + mu_ap(2);
det_Sigma = det(Sigma_ap); % Determinant of Sigma
inv_Sigma = inv(Sigma_ap); % Inverse of Sigma

for i = 1:N
    for j = 1:N
        x_vec = [X_shifted(i,j); Y_shifted(i,j)]; % Vector [x - mu_x, y - mu_y]
        rho_ap(i,j) = exp(-0.5 * x_vec' * inv_Sigma * x_vec) / (2*pi * sqrt(det_Sigma));
    end
end

figure; hold on
% surf(x_padded, x_padded, abs(C));
% xlabel('X'); ylabel('Y'); title('Convolution Result'); 
% colorbar;shading flat;
surf(X_original, Y_original, rho_ap);

%%

