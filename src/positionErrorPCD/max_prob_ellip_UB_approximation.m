function prob = max_prob_ellip_UB_approximation(s1, s2, Sigmax)
% Get bounding ellipsoid for the Minkowski sum of s1 and s2
[~, Sigmaf] = get_bounding_ellip(s1, s2);

% In this setting, we let mx = m2-m1 and mf = 0
mx = s2.tc-s1.tc;
% mf = zeros(3,1);

% eigen value decomposition of Sigmax
% Decompose Sigmax to a diagonal matrix  lamda(independent variable)
[R, lamda] = eig(Sigmax);
% standard deviation
sigma1 = sqrt(lamda(1,1));
sigma2 = sqrt(lamda(2,2));
sigma3 = sqrt(lamda(3,3));

% space transformation R'
mx_transformed = R' * mx;
% new mean value
mu1 = mx_transformed(1);
mu2 = mx_transformed(2);
mu3 = mx_transformed(3);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!This may be wrong, the correct one should be
% R' * Sigmaf * R. 
Sigmaf_transformed = R * Sigmaf * R';

% Set sampled points numer n1 n2 n3 = 10
points_number = 10;
gauss_hermit_poly_roots = roots(hermite(points_number));

prob = 0;

for j1=1:points_number
    for j2=1:points_number
        for j3=1:points_number
            z1_j1 = gauss_hermit_poly_roots(j1);
            z2_j2 = gauss_hermit_poly_roots(j2);
            z3_j3 = gauss_hermit_poly_roots(j3);
            
            x = [sqrt(2)*sigma1*z1_j1+mu1; sqrt(2)*sigma2*z2_j2+mu2; sqrt(2)*sigma3*z3_j3+mu3];
            
            prob = prob + gauss_hermite_weight(points_number, z1_j1) * ...
                gauss_hermite_weight(points_number, z2_j2) * ...
                gauss_hermite_weight(points_number, z3_j3) * indicator_outer_ellipsoid(x, Sigmaf_transformed);
        end
    end
end

prob = prob * pi^(-3/2);
    
end

function y = indicator_outer_ellipsoid(x, Q)
if x' / Q * x <=1 
   y=1;
else
    y=0;
end
end

function omega = gauss_hermite_weight(n, z_j)
omega = (2^(n-1)*factorial(n)*sqrt(pi)) / (n^2*hermiteH(n-1, z_j)^2);
end
