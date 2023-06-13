function x = get_vee_vector(G)
%GET_VEE_VECTOR Computer the vee vector of a 4by4 homogeneous transfermation
%matrix
%Inputs:
%   G: 4by4 homogeneous transfermation matrix of SE(3) group
%Outputs:
%   x: vee vector in se(3) space

% Get head vector of rotation R
theta = acos((trace(G(1:3,1:3))-1)/2);

if isreal(theta) == 0
    warning("Rotation angle is a complex number in function get_vee_vector()");
elseif theta==0
        warning("Rotation is an indentity matrix");
        R_head_vector = zeros(3,1);
        x = zeros(6,1);
        x(1:3,1) = R_head_vector;
        x(4:6,1) = G(1:3,4);
        return
end

R_skew_matrix = theta*(G(1:3, 1:3) - G(1:3,1:3)') / (2*sin(theta));
R_head_vector = vex(R_skew_matrix);

% Get corresponding translation part 
t = G(1:3,4);
v = left_jacob_inv(R_head_vector)*t;

x = zeros(6,1);
x(1:3,1) = R_head_vector;
x(4:6,1) = v;

end


function jl_inv = left_jacob_inv(x) 
    jl_inv  = eye(3) - 1*skew(x)/2 + (1/norm(x)^2 - (1+cos(norm(x)))/(2*norm(x)*sin(norm(x))))*skew(x)*skew(x);

end