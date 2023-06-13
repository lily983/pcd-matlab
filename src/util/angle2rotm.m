function R = angle2rotm(theta)
% Get the rotation matrix in 2D paramterized by the rotation angle

R = [cos(theta), -sin(theta);
    sin(theta), cos(theta)];
end