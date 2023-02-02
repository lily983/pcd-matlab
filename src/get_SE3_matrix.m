function G = get_SE3_matrix(vee)
%GET_SE3_matrix Get the matrix from the vee vector (exponential)
%Input
%   vee: vee vector of se(3)
%Output
%   G: 4by4 matrix of SE(3)

G = zeros(4,4);

x = vee(1:3,1);
G(1:3,1:3) = eye(3) + sin(norm(x))*skew(x)/norm(x) + (1-cos(norm(x)))*skew(x)*skew(x)/sqrt(norm(x));

G(1:3,4) = left_jacob(vee(1:3,1))*vee(4:6,1);
G(4,4) = 1;

end

function jl = left_jacob(x)
    jl = eye(3) + (1-cos(norm(x)))/(norm(x)^2)*skew(x) + (norm(x) - sin(norm(x)))/(norm(x)^3)*skew(x)*skew(x);
end
