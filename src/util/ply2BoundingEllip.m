function ellip = ply2BoundingEllip(ply, plotEllip)
if nargin == 1
    plotEllip = 0;
end

% Load the PLY file
ptCloud = pcread(ply);

% Downsample the point cloud using random sampling
sampleFraction = 0.2; % Adjust the fraction to keep 10% of the points
ptCloudDownsampled = pcdownsample(ptCloud, 'random', sampleFraction);

% Extract the downsampled 3D points
points = ptCloudDownsampled.Location;

% % Extract 3D points
% points = ptCloud.Location; % This gives an Nx3 matrix of 3D points

[invSigma1, e1] = MinVolEllipse(points', 0.001);
Sigma1 = inv(invSigma1);
[R, lamda_power2, ~] = svd(Sigma1);
axang_ellip = rotm2axang(R);

ellip(1,:) = [sqrt(diag(lamda_power2))', 1, 1, e1', axang_ellip(1, 1:3), axang_ellip(1, 4)];

if plotEllip
    % plot bounding ellips 
    figure; hold on
%         plot_ellipse(Sigma1, e1, 'edgecolor', 'g');
    sq_ellip = SuperQuadrics({sqrt(diag(lamda_power2))', [1, 1], [0,0], e1, rotm2quat(R), [20,20]});
    sq_ellip.PlotShape('y',0.5);
    scatter3(points(:,1), points(:,2), points(:,3), '.');
end

% get csv file for bounding ellipsoids
ellip_table = table(ellip);
writetable(ellip_table, 'bounding_ellip_for_ply.csv');

end