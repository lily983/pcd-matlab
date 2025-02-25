function A_new = closedFormSurfaceEllip(A, R, idea)
%idea: 1, 2, 3, 4
% A: axes, diagonal matrix
% R: rotation samples 3x3xN
% A_new: new encapsulating surface

N = size(R,3);

if idea == 1
    A_new = zeros(3,3);
    for i = 1:N
        A_new = A_new + R(:,:,i) * A;
    end
    A_new = A_new ./ N;
elseif idea == 2
    A_new = zeros(3,3);
    for i = 1:N
        A_new = A_new + R(:,:,i) * A * R(:,:,i)';
    end
    A_new = A_new ./ N;
elseif idea == 3
    
end