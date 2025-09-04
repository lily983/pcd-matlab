function visualize_pose_error(s1, s2)
N=[20, 20];


mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
mu_vee = get_vee_vector(mu_g);

Sigma_g =8.0000e-08*eye(6);
Sigma_g(4:6,4:6) =  8.0000e-03*eye(3);

n = 20;
gi_vee_rand = mvnrnd(mu_vee', Sigma_g, n);

s3 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});

for i=1:size(gi_vee_rand, 1)
    gi_matrix = get_SE3_matrix(gi_vee_rand(i,1:6)');
    gi_tc = gi_matrix(1:3,4);
    gi_q = rotm2quat(gi_matrix(1:3,1:3));
    s3.tc = gi_tc;
    s3.q = gi_q;
    s3.PlotShape(hex2rgb('45AC59'), 0.2);
    pause(0.5)
end