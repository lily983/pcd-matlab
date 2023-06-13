function [a1,a2,b1,b2] = get_gmm_param(Sigma, Ratio, mb)
% Check dimension 
dimension = size(Sigma,2);

b1 = mb;
b2 = Ratio*mb;

syms ma;
a1 = ma*Ratio;
a2 = ma*-1;
eqn = a1 * exp(-b1/2) + a2 * exp(-b2/2) == (2*pi)^(dimension/2)*sqrt(det(Sigma));
ma = double(solve(eqn, ma, 'Real', true));

a1 = ma*Ratio;
a2 = ma*-1;
end