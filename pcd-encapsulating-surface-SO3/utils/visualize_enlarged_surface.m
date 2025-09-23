function visualize_enlarged_surface(obj, R, x_new)

figure; hold on; axis equal

for i=1:size(R,3)
    obj_i = SuperQuadrics({obj.a, obj.eps, [0,0], zeros(3,1), rotm2quat(R(:,:,i)), [20,20]});
    obj_i.PlotShape('b', 0.05,0.3);
end
obj.PlotShape('y', 0.4,0.2);

plotSurface(x_new, obj.N, 'r', 0.4); % red
end