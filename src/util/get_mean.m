function mu = get_mean(g)
% Compute mean of a set of SE(3)
% param g: The set of SE(3) elements, each row is one element, column ranks
% from position to quaternion in order of x y z x y z w
% return: Mean of the set

% number of poses
num = size(g,1);

% Mean: initialization by  taking average of vector in se(3)
mu = eye(4);

mu_log = zeros(6,1);
for j=1:num
    T = zeros(4,4);
    T(1:3, 4) = g(j, 1:3)';
    T(1:3, 1:3) = quat2rotm([g(j, 7), g(j, 4:6)]);
    mu_log = mu_log + get_vee_vector(T);
end
mu = get_SE3_matrix(mu_log./num);

% Iteratively get mean
mu_log = ones(6,1);
max_num = 10;
tol = 1e-05;
count = 1;

while norm(mu_log)>=tol && count <=max_num
    mu_log = zeros(6,1);
    for i=1:num
        T = zeros(4,4);
        T(1:3, 4) = g(j, 1:3)';
        T(1:3, 1:3) = quat2rotm([g(j, 7), g(j, 4:6)]);
        d_g_log = get_vee_vector(mu \ T);
        mu_log = mu_log + d_g_log;
    end 
    mu = mu * get_SE3_matrix(mu_log./num);
    count = count +1;
    norm(mu_log)
end

% if count > max_num
%     error("failed to compute mean");
% end

end