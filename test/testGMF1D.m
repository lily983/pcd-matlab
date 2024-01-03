clc;clear all; close all

% n: generate n random Sigma
n = 100;
Sigma = rand(n,1) * 5;

% generate x samples
x = linspace(-15, 15, 1000);

% here we ignore the center 
mu = 0;

% record non-upper bound case
error_a = {};
error_b = {};
error_Sigma = {};

% GMF 5SG sets
testProgress = waitbar(0, 'Starting test GMF1D upper bound for 5SG sets');
for i = 1:n
    [a, b] = get_gmm_param(Sigma(i), '5SG');
    
    f = GMF(mu, Sigma(i), a, b, x');
    
    iota = rectangularPulse(-sqrt(Sigma(i)), sqrt(Sigma(i)), x);
    
    difference_5SG = f - iota;

    if all(difference_5SG>=0)==false
        error_a = [error_a; a];
        error_b = [error_b; b];
        error_Sigma = [error_Sigma; Sigma(i)];
        error('Find non-upper bound')
    end
    waitbar(i/n, testProgress, sprintf('Progress(GMF1D, 5SG): %d %%', floor(i/n*100)));
end
close(testProgress);

% GMF 2SG sets
testProgress = waitbar(0, 'Starting test GMF1D upper bound for 2SG sets');
for i = 1:n
    [a, b] = get_gmm_param(Sigma(i), '2SG');
    
    f = GMF(mu, Sigma(i), a, b, x');
    
    iota = rectangularPulse(-sqrt(Sigma(i)), sqrt(Sigma(i)), x);
    
    difference_2SG = f - iota;

    if all(difference_2SG>=0)==false
        error_a = [error_a; a];
        error_b = [error_b; b];
        error_Sigma = [error_Sigma; Sigma(i)];
        error('Find non-upper bound')
    end
    waitbar(i/n, testProgress, sprintf('Progress(GMF1D, 2SG): %d %%', floor(i/n*100)));
end
close(testProgress);
%% Function
function y = GMF(mu, Sigma, a, b, x)
sizeParameters = max(size(a));

y = zeros(1, max(size(x)));

for i = 1:sizeParameters
    for j = 1:max(size(x))
        y(j) = y(j) + a(i) * b(i)^(1/2) * exp( (-b(i)/2) * (x(j)-mu)' / Sigma * (x(j)-mu));
    end
end

y = y ./ ( (2*pi)^(1/2) * det(Sigma)^0.5);

end
