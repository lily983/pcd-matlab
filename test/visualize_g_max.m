function visualize_g_max(s1,s2,g_max)

s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

s3.tc = g_max(1:3,4);
s3.q = rotm2quat(g_max(1:3,1:3));

figure; hold on;
s1.PlotShape('b',0.2);
s2.PlotShape('g',0.2);
s3.PlotShape('r',0.1);


end