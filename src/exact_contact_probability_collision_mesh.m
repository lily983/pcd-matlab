function prob = exact_contact_probability_collision_mesh(mu, Sigma, s1, s2, N)


prob = 0;

mu_vee = get_vee_vector(mu);

samples = mvnrnd(mu_vee', Sigma, N);

s3 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});
% 
% figure; hold on;
% s1.PlotShape('b', 0.3);
% s2.PlotShape('g', 0.1);

for i=1:N
    i
    g_rand = get_SE3_matrix(samples(i, 1:6)');
    
    % Update pose error
    s3.tc = g_rand(1:3,4);
    s3.q = rotm2quat(g_rand(1:3,1:3));
    
   % for two sphere objects, use simplier method to test if they collide
    if s1.eps(1)==1 && s1.eps(2)==1 && s3.eps(1) ==1 && s3.eps(2)==1 && isequal(s1.a./s1.a(1), ones(1,3))==true && isequal(s3.a./s3.a(1), ones(1,3))==true
       if norm(s1.tc-s3.tc)<=(s1.a(1)+s3.a(1))
           prob = prob+1;      
       end
    else
       if collision_mesh(s1,s3)
           prob = prob+1;
%            s3.PlotShape('r',0.1)
%            pause(0.1)
       end
   end
end
prob=prob/N;
end