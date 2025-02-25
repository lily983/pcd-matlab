%% Loading rotations from samples
R = csv2group(samples, 'SO3');

%  R_filter = filterData(R, mu_SO3);
R_filter = R;

%% Loading obj parameters
axes = [0.034,0.021, 0.016];
obj = SuperQuadrics({axes, [0.2,0.2], [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

%% Get encapsulating surface points
% x_new_1 = idea1(obj, R_filter);
% 
% x_new_2 = idea2(obj, R_filter);
% 
% x_new_3 = idea3(obj, R_filter);

x_new_4_13 = idea4scale(obj, R_filter, size(R_filter, 3)/1.3);

x_new_4_11 = idea4scale(obj, R_filter, size(R_filter, 3)/1.1);

x_new_8_11 = idea8(obj, R_filter, size(R_filter, 3)/1.1);

x_new_8_13 = idea8(obj, R_filter, size(R_filter, 3)/1.3);
% [x_new_51, d_51] = idea5(obj, R_filter, '1' );
% 
% [x_new_52, d_52] = idea5(obj, R_filter, '2' );

%% 
figure; hold on; axis equal
plotSurface(x_new_4_13, obj.N,  hex2rgb('45AC59'), 0.1); % dark green
plotSurface(x_new_4_11, obj.N, hex2rgb('98646b'), 0.1);% grape
plotSurface(x_new_8_11, obj.N, hex2rgb('2023c7'), 0.1); % bright purple
% plotSurface(x_new_4, obj.N, hex2rgb('EBBF00')); % ginger yellow
plotSurface(x_new_8_13, obj.N, hex2rgb('ff0010'), 0.1); % red
% plotSurface(mu_SO3*x_new_52, obj.N, hex2rgb('45AC59')); %  dark green


%% Plot obj at mean rotation and its slightly rotated copies 
close all;
figure; hold on; axis equal
% obj_mean = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});
% obj_mean.PlotShape('b', 0.1,0.1);

[R_extreme, d_R] = getExtremeR(R_filter, mu_SO3, Sigma_SO3);

for i=1:size(R_extreme,3)
%     obj_i = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(R(:,:,randi([1, size(R,3)]))), [20,20]});
    obj_i = SuperQuadrics({axes, [0.2,0.2], [0,0], zeros(3,1), rotm2quat(R_extreme(:,:,i)), [20,20]});
    obj_i.PlotShape('b', 0.001,0.15);
    pause(0.1);
end
