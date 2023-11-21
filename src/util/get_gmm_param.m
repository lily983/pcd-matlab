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
    case '1SG'
        a= [2e+09];
        b= [40.9949];
    case '5SG'
        a=[500 500 500 500 500];
        b=[12.2474, 14.7835, 13.0516, 18.6186, 22.1856];
end

% Get the coeffecient ma
ma=0;

for i=1:size(a, 2)
    ma=ma+a(i)*exp(-b(i)/2);
end

ma = ((2*pi)^(dimension/2)*sqrt(det(Sigma)))/ma;
a = ma.*a;

end