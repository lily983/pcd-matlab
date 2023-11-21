function prob = max_prob_ellip_UB_exact(s1, s2, Sigmax)
% Check dimension 
% dimension = size(Sigmax, 2);

% Get bounding ellipsoid for the Minkowski sum of s1 and s2
[~, Sigmaf] = get_bounding_ellip(s1, s2);

% In this setting, we let mx = m2-m1 and mf = 0
mx = s2.tc-s1.tc;
% mf = zeros(3,1);

% eigen value decomposition of Sigmaz
[R, lamda] = eig(Sigmax^0.5 * Sigmaf * Sigmax^0.5);

% Get variable b
b=R' / Sigmax^0.5 * mx;

prob = fixed_point(b, diag(lamda));
if imag(prob)>1e-02
    prob = NaN;
else
    prob = real(prob);
end

end


function prob = fixed_point(b, lamda)
prob=0;

% matlab array starts from position 1 instead of 0. So k starts from 1
% while in paper it starts from position 0
k=1;
c0 = exp(-0.5 * (b(1)^2 + b(2)^2 + b(3)^2)) * ((2*lamda(1)) * (2*lamda(2)) * (2*lamda(3)))^(-0.5);
c_list(k) = c0;

delta = (-1)^(k-1)*c_list(k)/gamma(3/2+k);
prob = prob+delta;

k=2;
c_next = fixed_point_iteration(c_list, b, k, lamda);

% keep getting more series until the increment is less than threshold
while k<30
    c_list(k) = c_next; 
    delta = (-1)^(k-1)*c_list(k)/gamma(3/2+k);
    prob = prob+delta;
    
    k=k+1;
    c_next = fixed_point_iteration(c_list, b, k, lamda);
end

end

function c_k = fixed_point_iteration(c, b, k, lamda)
c_k=0;

for i=0:1:(k-2)
    c_k = c_k + get_dk(b, k-i, lamda) * c(i+1);
end

c_k = (1.0/(k-1)) * c_k;
end

function d_k = get_dk(b, k, lamda)
d_k = 0;

for i=1:3
    d_k = d_k + (1 - (k-1)*b(i)^2)*(2*lamda(i))^(-(k-1));
end

d_k = d_k * 0.5;
end