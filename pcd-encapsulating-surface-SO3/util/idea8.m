function x_new = idea8(obj, R, scale)
% obj: superquadric
% R: rotation samples 3x3xN
% x_new: new encapsulating surface points 3 x m

% x_original  = obj.GetPoints();

x_new = zeros(3, obj.N(1) * obj.N(2));

% R_new = zeros(3,3);

n = obj.GetNormals();

% figure; hold on
% obj.PlotShape('b', 0.1,0.1);
% 
n1_new = R(:,:,1) * n;

for i = 1:size(R,3)
% x_new = x_new + R(:,:,i) * obj.GetPointsFromNormal(R(:,:,i)' * n);

x_new = x_new + R(:,:,i) * obj.GetPointsFromNormal(R(:,:,i)' * n1_new);

% plotSurface(R(:,:,i) * obj.GetPointsFromNormal( R(:,:,i)' * n1_new), obj.N, hex2rgb('45AC59'), 0.1);
% plotSurface(x_new, obj.N, hex2rgb('45AC59'), 0.1);
end

x_new = x_new  ./ scale;
    
end