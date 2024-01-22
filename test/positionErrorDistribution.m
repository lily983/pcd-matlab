function Sigma = positionErrorDistribution(distribution)

% if nargin==1
%     distribution='large';
% end

Sigma = zeros(3);
Sigma(1,1) = 7.0000e-00;
Sigma(2,2) = 8.0000e-00;
Sigma(3,3) = 9.0000e-00;

switch distribution
    case 'small'
        scaleValue = 4e-04;
    case 'large'
        scaleValue = 6e-03;
    case 'none'
        scaleValue = 0;
end

% if strcmp(objectType, 'sphere')
%     sizeScale = 0.1*sizeScale;
% end

Sigma = scaleValue.*Sigma;

end