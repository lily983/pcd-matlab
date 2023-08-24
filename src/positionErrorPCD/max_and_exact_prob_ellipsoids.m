function prob = max_and_exact_prob_ellipsoids(Sigma, s1, s2, method)
% Get the maximum and exact value of PCD for ellipsoids with Gaussian
% position uncertainty (benchmark method)
% 'method': 'max' or 'exact'

prob=0;

% Difference between centers of two ellipsoids
mu_y = s2.tc - s1.tc;
Sigma_y = Sigma;

% Defining matrices
B = quat2rotm(s1.q) * diag(s1.a).^2 * quat2rotm(s1.q)';
C = quat2rotm(s2.q) * diag(s2.a).^2 * quat2rotm(s2.q)';

% Matrix M
C_bar = sqrtm(B) / C * sqrtm(B);
c_bar = sqrtm(B) \ (s2.tc - s1.tc);

M=zeros(6,6);
M(1:3, 1:3) = inv(C_bar);
M(1:3, 4:6) = -1*eye(3);
M(4:6, 1:3) = -1*(sqrtm(C_bar) \ c_bar)*(sqrtm(C_bar) \ c_bar)';
M(4:6, 4:6) = inv(C_bar);

% Get the minimal eigen value of M
lamda_0 = min(real(eig(M)));

% Matrix A
D = sqrtm(B) / (lamda_0*eye(3) - inv(C_bar)) / sqrtm(B);
A = D' / B * D;
 
switch method
    case 'exact'
        transform_A = sqrtm(Sigma_y) * A * sqrtm(Sigma_y);
        [P, lamda] = eig(transform_A);
        vec_b = P' / sqrtm(Sigma_y) * mu_y;
        if norm(eig(P*P') - ones(3,1))>1e-02
            prob = NaN;
            return;
        end
        prob = fixed_point(lamda_0, vec_b, diag(lamda));
        if imag(prob)>1e-02
            prob = NaN;
        else
            prob = real(prob);
        end
    case 'max'
        %compute beta
        [expectation, standard_deviation] = from_quadratic(mu_y, Sigma_y, A);
        beta = expectation + standard_deviation;
        if beta < 1/lamda_0^2
            prob=1;
        else
            prob = (beta-expectation) / (beta-1/lamda_0^2);
        end
end
end

function prob = fixed_point(lamda_0, b, lamda)
prob=0;
v = 1/lamda_0^2;

k=1;
c0 = exp(-0.5 * (b(1)^2 + b(2)^2 + b(3)^2)) * ((2*lamda(1)) * (2*lamda(2)) * (2*lamda(3)))^(-0.5);
c_list(k) = c0;

delta = (-1)^(k-1)*c_list(k)*v^(3/2+k-1)/gamma(3/2+k-1+1);
prob = prob+delta;

k=2;
c_next = fixed_point_iteration(c_list, b, k, lamda);

while abs(delta)>1e-03
    c_list(k) = c_next; 
    delta = (-1)^(k-1)*c_list(k)*v^(3/2+k-1)/gamma(3/2+k-1+1);
    prob = prob+delta;
    
    k=k+1;
    c_next = fixed_point_iteration(c_list, b, k, lamda);
end

end

function c_k = fixed_point_iteration(c, b, k, lamda)
c_k=0;
for i=0:1:(k-2)
    c_k = c_k + (1.0/(k-1)) * get_dk(b, k-1-i, lamda) * c(i+1);
end
end

function d_k = get_dk(b, k, lamda)
d_k = 0;

for i=1:3
    d_k = d_k + (1 - k*b(i)^2)*(2*lamda(i))^(-k);
end

d_k = d_k * 0.5;
end

function [expectation, standard_deviation] = from_quadratic(mu, Sigma, A)
% sample x from its distribution
N=1e+04;
samples = mvnrnd(mu, Sigma, N);
quadratic_form = zeros(1, N);
i=1;
while abs(mean(quadratic_form) - trace(A*Sigma)-mu'*A*mu)>1e-2 && i<N
    i=i+1;
    quadratic_form(i) = samples(i,:) * A * samples(i,:)';
end
expectation = trace(A*Sigma)+mu'*A*mu;
standard_deviation = std(quadratic_form);
end

