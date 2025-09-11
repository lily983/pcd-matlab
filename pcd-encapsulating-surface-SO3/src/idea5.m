function [x_new, d, d_index] = idea5(obj, R, methodName)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m
A = diag(obj.a);

[a, index] = max(obj.a);

if index==1
    ei = [1;0;0];
elseif index ==2
    ei = [0;1;0];
else
    ei = [0;0;1];
end

results = zeros(1, size(R,3));

for i = 1:size(R,3)
    if strcmp(methodName, '1')
        results(i) = a - 1/norm(A\R(:,:,i)*ei);
    elseif strcmp(methodName, '2')
%         results(i) = a - 1/norm(A\R(:,:,i)*ei);
         results(i) = (a - 1/norm(A\R(:,:,i)*ei)) * ( ei' * R(:,:,i)' / A * R(:,:,i) * ei ) / norm(A\R(:,:,i)*ei);
%     elseif strcmp(methodName, '3')
%         results(i) = (ei' * (R(:,:,i)' / A * R(:,:,i) *A - (R(:,:,i)' / A * R(:,:,i)) ./ norm(A\R(:,:,i)*ei) ) * ei) ./ norm(A\R(:,:,i)*ei) ;
%         results(i) = (a - 1/norm(A\R(:,:,i)*ei)) * ( ei' * R(:,:,i)' / A * R(:,:,i) * ei ) / norm(A\R(:,:,i)*ei);
    end
end

[d, d_index] = max(results);

% x_original  = obj.GetPoints();

n = obj.GetNormals();

% x_new = zeros(3, obj.N(1) * obj.N(2));

x_new = (A^2 *n) ./ vecnorm(A * n, 2, 1) + d .* n;

% x_new = x_original + (d .* A^2 *n) ./ vecnorm(A^2 * n, 2, 1);
    
end