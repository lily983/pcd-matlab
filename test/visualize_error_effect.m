close all; clear; clc;
add_path()

N=[20, 20];
s1 = SuperQuadrics({0.02*ones(1,3), [0.8,0.7], [0, 0]...
zeros(3,1), [1, 0, 0, 0], N});

init_quat = rand(1,4);
init_quat = init_quat/norm(init_quat);
init_posi = 0.08*rand(3,1);

s2 = SuperQuadrics({[0.01,0.02,0.03], [0.4,0.8], [0, 0]...
        init_posi, init_quat, N});
%%
figure; hold on; axis equal;
lightangle(gca,45,30);
lighting gouraud;

s1.PlotShape('g', 0.8);
s2.PlotShape('y', 0.8);
mu_g =  [quat2rotm(init_quat), init_posi; 0, 0, 0, 1];
mu_vee = get_vee_vector(mu_g);

Sigma_g =10.0000e-06*eye(6);
Sigma_g(4:6,4:6) =  5.0000e-06*eye(3);
[prob_max, g_max] = max_contact_probability(mu_g, Sigma_g, s1, s2)
s3 = SuperQuadrics({[0.01,0.02,0.03], [0.4,0.8], [0, 0]...
        g_max(1:3,4), rotm2quat(g_max(1:3, 1:3)), N});
s3.PlotShape('r',0.4)
%%
m1 = s1.GetGradientsCanonical();
mink = MinkSumClosedForm(s1,s3,quat2rotm(s1.q),quat2rotm(s3.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));

surf(X, Y, Z,...
 'EdgeColor', 'b', 'EdgeAlpha', 0.3,...
 'FaceAlpha', 0.1, 'FaceColor', 'b'); %plot contact space

%%
mu_g =  [quat2rotm(init_quat), init_posi; 0, 0, 0, 1];
mu_vee = get_vee_vector(mu_g);

Sigma_g =10.0000e-06*eye(6);
Sigma_g(4:6,4:6) =  5.0000e-06*eye(3);

n = 100;
gi_vee_rand = mvnrnd(mu_vee', Sigma_g, n);

for i=1:size(gi_vee_rand, 1)
    gi_matrix = get_SE3_matrix(gi_vee_rand(i,1:6)');
    gi_tc = gi_matrix(1:3,4);
    gi_q = rotm2quat(gi_matrix(1:3,1:3));
    s2.tc = gi_tc;
    s2.q = gi_q;
    s2.PlotShape('y', 0.1);
end
