function g = csv2group(csv, groupName)
% csv: each row is x y z x y z w
% groupName: SO3 PCG3
% g: 4 x 4 x N homogeneous transformation matrices 
% g: 3 x 3 x N rotation matrices

N = size(csv, 1);

if nargin == 1
    groupName = 'PCG3';
end

if strcmp(groupName, 'PCG3')
    g = zeros(4, 4, N);
    for i = 1:N
        g(1:3,4,i) = csv(i, 1:3);
        g(1:3, 1:3, i) = quat2rotm([csv(i, 7), csv(i, 4:6)]);
        g(4, 4, i)=1;
    end
elseif strcmp(groupName, 'SO3')
    g = zeros(3, 3, N);
    for i = 1:N
        g(1:3, 1:3, i) = quat2rotm([csv(i, 7), csv(i, 4:6)]);
    end
elseif strcmp(groupName, 'R3')
    g = zeros(4, 4, N);
    for i = 1:N
        g(1:3,4,i) = csv(i, 1:3);
        g(1:3, 1:3, i) = eye(3);
        g(4, 4, i)=1;
    end
end
    

end