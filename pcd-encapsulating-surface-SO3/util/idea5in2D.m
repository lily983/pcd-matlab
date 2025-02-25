function [x_new, d, d_index, d_max, theta_max] = idea5in2D(obj, R, methodName)
% obj: SuperEllipse
% R: rotation samples 2x2xN
% x_new: new encapsulating surface points 2 x m
A = diag(obj.a);
e1 = [1;0];
a  = obj.a(1);

results = zeros(1, size(R,2));

for i = 1:size(R,3)
    if strcmp(methodName, '1')
        results(i) = a - 1/norm(A\R(:,:,i)*e1);
%          results(i) = (e1' * R(:,:,i)' / A * R(:,:,i) * A * e1 - e1' * R(:,:,i)' / A * R(:,:,i) * e1/norm(A \ R(:,:,i) * e1))  / norm(A \ R(:,:,i) * e1);
%         results_2(i) = (a - 1/norm(A \ R(:,:,i) * e1)) * ( e1' * R(:,:,i)' / A^2 * R(:,:,i) * e1 ) / norm(A^2 \ R(:,:,i) * e1);
    elseif strcmp(methodName, '2')
        % projection of d in the direction of n
        results(i) = (e1' * R(:,:,i)' / A * R(:,:,i) * A * e1 - e1' * R(:,:,i)' / A * R(:,:,i) * e1/norm(A \ R(:,:,i) * e1))  / norm(A \ R(:,:,i) * e1);
    elseif strcmp(methodName, '3')
%         results(i) = (a - 1/norm(A \ R(:,:,i) * e1)) * ( e1' * R(:,:,i)' / A^2 * R(:,:,i) * e1 ) / norm(A^2 \ R(:,:,i) * e1);
        results(i) = a - 1/norm(A\R(:,:,i)*e1);
    end
end

[d, d_index] = max(results);

m = obj.GetGradients();

n = m ./ vecnorm(m, 2, 1);

if strcmp(methodName, '3')
    theta = -1*(0:0.01:pi/2);

    for j=1:size(theta,2)
        R_theta = rot2(theta(j));
        e_theta = R_theta * e1;
        d_theta(j) = abs(1/norm(R(:,:,d_index)/A*R(:,:,d_index)'*e_theta) - 1/norm(A \ e_theta));
    end

    [d_max , theta_index] = max(d_theta);
    theta_max = theta(theta_index);
    
    x_new = obj.GetPoints() + d_max .* n;
    return
end

d_max = nan;
theta_max = nan;
x_new = obj.GetPoints() + d .* n;
   
end