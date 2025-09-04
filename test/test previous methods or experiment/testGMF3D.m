clc;clear all; close all

% n: generate n random Sigma
n = 100;
semi_axes_square = rand(3, n)*5;
quatEllip = randrot(1, n);

% generate samples
% Because mu=0 and all functions are symmetric about xz, yz-planes, only
% computes points from positive x, y-axes.
x = linspace(0, 15, 1000);
y = linspace(0, 15, 1000);
z = linspace(-15, 15, 2000);

% here we ignore the center 
mu = zeros(3,1);

% record non-upper bound case
error_a = {};
error_b = {};
error_Sigma = {};

% GMF 5SG sets
testProgress = waitbar(0, 'Starting test GMF3D upper bound for 5SG sets');
for i = 1:n
    Sigma_matrix = quat2rotm(quatEllip(1, i))' * diag(semi_axes_square(:,i)) * quat2rotm(quatEllip(1,i));
    
    [a, b] = get_gmm_param(Sigma_matrix, '5SG');
    
    f = GMF3D(mu, Sigma_matrix, a, b, x, y, z);
    
    iota = indicatorFunction(mu, Sigma_matrix, x, y, z);
    
    difference = matrixDiff(iota, f);
    
    if any(difference<0)
        error_a = [error_a; a];
        error_b = [error_b; b];
        error_Sigma = [error_Sigma; Sigma(i)];
        error('Find non-upper bound')  
    end
    
    waitbar(i/n, testProgress, sprintf('Progress(GMF3D, 5SG): %d %%', floor(i/n*100)));
end
close(testProgress);

% GMF 2SG sets
testProgress = waitbar(0, 'Starting test GMF3D upper bound for 2SG sets');
for i = 1:n
    Sigma_matrix = quat2rotm(quatEllip(1, i))' * diag(semi_axes_square(:,i)) * quat2rotm(quatEllip(1,i));
    
    [a, b] = get_gmm_param(Sigma_matrix, '2SG');
    
    f = GMF3D(mu, Sigma_matrix, a, b, x, y, z);
    
    iota = indicatorFunction(mu, Sigma_matrix, x, y, z);
    
    difference = matrixDiff(iota, f);
    
    if any(difference<0)
        error_a = [error_a; a];
        error_b = [error_b; b];
        error_Sigma = [error_Sigma; Sigma(i)];
        error('Find non-upper bound')  
    end
    
    waitbar(i/n, testProgress, sprintf('Progress(GMF3D, 2SG): %d %%', floor(i/n*100)));
end
close(testProgress);
%% Functions
% Indicator function in 3D: for ellipsoidal point sets
function result = indicatorFunction(mu, Sigma, x, y, z)
% If points (x;y;z) is inside the ellipsoid, (x-mu)' * Sigma^-1 * (x-mu) >=1, then
% function returns 1
sizePoints = max(size(x));

result = ones(sizePoints, sizePoints, sizePoints);

for i=1:sizePoints
    for j=1:sizePoints
        for k=1:sizePoints
            point = [x(i); y(j); z(k)];
            if (point-mu)' / Sigma * (point-mu) > 1
                result(i, j, k) = 0;
            end
        end
    end
end

end

function result = GMF3D(mu, Sigma, a, b, x, y, z)
% Get the function value at points (x;y;z) of function GMF
sizePoints = max(size(x));

result = zeros(sizePoints, sizePoints, sizePoints);

sizeParameters = max(size(a));

for i=1:sizePoints
    for j=1:sizePoints
        for k=1:sizePoints
            point = [x(i); y(j); z(k)];
            for m=1:sizeParameters
                result(i,j,k) = result(i,j,k) + a(m) * b(m)^(3/2) * exp( (-b(m)/2) * (point-mu)' / Sigma * (point-mu));
            end
            result(i,j,k) = result(i,j,k) ./ ( (2*pi)^(3/2) * det(Sigma)^0.5);
        end
    end
end

end

% Get the matrix difference Y-X along z axis 
function result = matrixDiff(X, Y)
% Size of matrices X and Y is the same
result = zeros(size(X,1) * size(X,2) * size(X,3), 1);
for i=1:size(X,1)
    for j=1:size(X,2)
        for k=1:size(X,3)
           result(i+j+k) =  Y(i,j,k)-X(i,j,k);
        end
    end
end
end

