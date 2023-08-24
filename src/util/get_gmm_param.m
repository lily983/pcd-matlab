function [a, b] = get_gmm_param(Sigma, method)
%%Get the parameters for the Gaussian mixture function
% Sigma: the defining matrix of the ellipsoid
% method: which set of parameters to use
% a, b: parameters of GMM function

% Check dimension 
dimension = size(Sigma,2);

% Check inputs. Set default value
if nargin==1
    method='5SG';
end

% Choose methods: number of single gaussians
switch method
    case '2SG'
        a=[2.9436   15.4697];
        b=[30.0000    3.6882];
    case '3SG'
        a=[-9.5616 7.9123 4.4867];
        b=[7.1483 3.4981 2.6616];
    case '4SG'
        a=[-10   -10    20    20];
        b=[13.4603   13.4603    6.6984    4.8889];
    case '5SG'
        a=[3.9302   40.0000   -1.3492  -19.9206    7.5397];
        b=[22.6667    6.0476    4.3333   30.0000    5.1428];
end

% Get the coeffecient ma
ma=0;

for i=1:size(a, 2)
    ma=ma+a(i)*exp(-b(i)/2);
end

ma = ((2*pi)^(dimension/2)*sqrt(det(Sigma)))/ma;
a = ma.*a;

end