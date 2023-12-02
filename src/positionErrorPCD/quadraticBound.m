function [prob, t] = quadraticBound(s1, s2, mx, Sigma)
% This function reproduce methods using Markov inequality of random
% variable to get PCD value. The code is written based on 
% Eq. 22 in paper "Exact and Bounded Collision Probability for Motion
% Planning under Gaussian Uncertainty" 2022 RAL Antony Thomas
    % Inputs:
        % s1, s2: sphere or ellipsoid
        % Sigma: covariance matrix of position error distribution
    % Outputs:
        % prob: PCD value
        % t: computation time
tic;
prob=0;

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
 
%compute beta
[expectation, standard_deviation] = from_quadratic(mx, Sigma, A);
beta = expectation + standard_deviation;
if beta < 1/lamda_0^2
    prob=1;
else
    prob = (beta-expectation) / (beta-1/lamda_0^2);
end

t=toc;
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

