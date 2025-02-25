function Sigma = get_covariance(g, mu)
% Get covariance of a set of SE(3)
% param g: The set of SE(3) elements
% return: Covariance of the set
if nargin==1
    mu = get_mean(g);
end

num = size(g,1);

Sigma = zeros(6,6);
for i=1:num
    T = zeros(4,4);
    T(1:3, 4) = g(i, 1:3)';
    T(1:3, 1:3) = quat2rotm([g(i, 7), g(i, 4:6)]);
    y_i = get_vee_vector(mu \ T);
    Sigma = Sigma + y_i * y_i';
end
Sigma = Sigma ./ num;

end