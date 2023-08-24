function Sigma = positionErrorDistribution(objectType, distribution)

if nargin==1
    distribution='sparse';
end

Sigma = zeros(3);
Sigma(1,1) = 7.0000e-00;
Sigma(2,2) = 8.0000e-00;
Sigma(3,3) = 9.0000e-00;

switch distribution
    case 'concentrate'
        sizeScale = 8e-04;
    case 'sparse'
        sizeScale = 3e-03;
end

if strcmp(objectType, 'sphere')
    sizeScale = 0.1*sizeScale;
end

Sigma = sizeScale.*Sigma;

end