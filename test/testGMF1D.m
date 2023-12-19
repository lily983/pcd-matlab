clc;clear all; close all

% n: generate n random Sigma
n = 10;
Sigma = rand(n,1) * 5;

% generate x samples
x = -15:0.01:15;

% here we ignore the center 
mu = 0;

% record non-upper bound case
error_a = {};
error_b = {};
error_Sigma = {};

% GMF 5SG sets
for i = 1:n
    [a, b] = get_gmm_param(Sigma(i), '5SG');
    
    f = GMF(x', mu, Sigma(i), a, b, 1);
    
    iota = rectangularPulse(-sqrt(Sigma(i)), sqrt(Sigma(i)), x);
    
    difference_5SG = f - iota;

    if all(difference_5SG>=0)==false
        error_a = [error_a; a];
        error_b = [error_b; b];
        error_Sigma = [error_Sigma; Sigma(i)];
        error('Find non-upper bound')
    end
end

% GMF 2SG sets
for i = 1:n
    [a, b] = get_gmm_param(Sigma(i), '2SG');
    
    f = GMF(x', mu, Sigma(i), a, b, 1);
    
    iota = rectangularPulse(-sqrt(Sigma(i)), sqrt(Sigma(i)), x);
    
    difference_2SG = f - iota;

    if all(difference_2SG>=0)==false
        error_a = [error_a; a];
        error_b = [error_b; b];
        error_Sigma = [error_Sigma; Sigma(i)];
        error('Find non-upper bound')
    end
end


