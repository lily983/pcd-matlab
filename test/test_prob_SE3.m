%% Test 1: pure traslation case. prob_x_max, prob_x_s1, prob_x_polygon, exact_prob_x
close all;clear;clc;
add_path()

N = [20,20];
figure; hold on; axis equal;
lightangle(gca,45,30);
lighting gouraud;

for i=1:50
    
    s1 = SuperQuadrics({0.5*ones(1,3), [1,1], [0, 0]...
    zeros(3,1), [1, 0, 0, 0], N});

    s2 = SuperQuadrics({0.5*ones(1,3), [1,1], [0, 0]...
        [1.01, 0, 0], [1, 0, 0, 0], N});

    [flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
    
    mu=s2.tc;
    
    Sigma = zeros(3);
    Sigma(1,1) = 0.0001;
    Sigma(2,2) = 0.0009;
    Sigma(3,3) = 0.0025;
    
    %Sigma = 0.001*eye(3);
%{
    mu_g =  [quat2rotm(s2.q), s2.tc; 0, 0, 0, 1];
    Sigma_g =1*eye(6);
    Sigma_g(4:6,4:6) = 0.1*eye(3);
    %}
   [pdf_max_x, x_max] = max_probability_pure_translation(mu, Sigma, s1, s2);
    
   s2_copy = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});
  % [pdf_max_g, g_max] = max_contact_probability(mu_g, Sigma_g, s1, s2_copy, false);
   
   v =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
    prob_max_x = v * pdf_max_x;

 %   pkf_value = pkf_3d_sphere(s1, s2);
   % prob_max_g = pkf_value*pdf_max_g;
    
    prob_max_x1 = mvnpdf(s1.tc, mu, Sigma) * v;
%{
    g1 = eye(4,4);
    g1(1:3,1:3) = quat2rotm(s1.q);
    g1(1:3,4) = s1.tc;
    g1_matrix = logm(g1);
    g1_vector = [vex(g1_matrix(1:3,1:3)); g1_matrix(1:3,4)];
    
    g_max = mu_g;
    g_skew_matrix = logm(g_max);
    g_skew_vector = [vex(g_skew_matrix(1:3,1:3)); g_skew_matrix(1:3,4)];
    
    prob_g1 = mvnpdf(g1_vector', g_skew_vector', Sigma_g) * pkf_value;
   %}
   % prob_max_cs = sqrt(pkf_value)*sqrt(sqrt(2)/sqrt((2*pi)^3*abs(det(Sigma))));
 

    exact_prob_x = exact_contact_probability_pure_translation(mu, Sigma, s1, s2, 1e+05);
   % exact_prob_g = exact_contact_probability(mu_g, Sigma_g, s1, s2, 10000);
    
    prob_max_poly = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2);
    
     s1.PlotShape('b', 0.5);
     s2.PlotShape('g', 0.5);

     %{
        result(1,i) = double(flag);
        result(2,i) = dist;
        %{
        result(3,i) = pdf_max_x;
        result(4,i) = pdf_max_g;
        result(5,i) = pkf_value;
        result(6,i) = v;
        result(7,i) = prob_max_x;
        result(8,i) = prob_max_g;
        result(9,i) = prob_max_x1;
        result(10,i) = prob_max_cs;
        %}
        result(3,i) = exact_prob_x;
        %result(12,i) = exact_prob_g;
        %result(13,i) = prob_g1;
        result(4,i) = prob_max_poly;
     %}   
     
        result(1,i) = double(flag);
        result(2,i) = dist;
        result(3,i) = pdf_max_x;
        result(4,i) = prob_max_x;
        result(5,i) = prob_max_x1;
        result(6,i) = exact_prob_x;
        result(7,i) = prob_max_poly;
    
end

%% Store results
pathname = '~/prob-collision-matlab/test/result';
test = "compare_max_x_and_max_poly" + "_exclude_collide";
filename = test + ".mat";
fig = test +".fig";
save(fullfile(pathname, filename), 'result')
savefig(fullfile(pathname, fig))

%% Plot results
figure; hold on;

plot(result(7,:),'b')
%plot(result(4,:),'y')
plot(result(6,:),'g')
%%
plot(result(11,:)./result(8,:))

%% Box plot 
figure;
difference  = result(4,:)' - result(3,:)';
boxplot(difference)

%%
points_1 = s1.GetPointsCanonical()';
shp_1 = alphaShape(points_1);
tri_1 = alphaTriangulation(shp_1);

points_2 = s2.GetPointsCanonical()';
shp_2 = alphaShape(points_2);

mink_points = zeros(size(points_1,1)+size(points_2,1),3);
k=1;
for i =1:size(points_1,1)
    for j=1:size(points_2,1)
      mink_points(k,:) = points_1(i,:)+points_2(j,:);
      k=k+1;
    end
end
shp_mink = alphaShape(mink_points);

%%
polygon_area = 0;
for i=1:1:size(tri_1,1)
    v1=points_1(tri_1(i,1),:);
    v2=points_1(tri_1(i,2),:);
    v3=points_1(tri_1(i,3),:);
    
    polygon_area = polygon_area+1/2*norm(cross(v1-v2, v1-v3));
end
%%
m1 = s1.GetGradientsCanonical();
[flag, dist, pt_cls, condition] = collision_cfc(s1,s2,'constrained');
mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;

X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));
%%
figure; hold on;

surf(X, Y, Z,...
 'EdgeColor', 'k', 'EdgeAlpha', 0.2,...
 'FaceAlpha', 0.1, 'FaceColor', 'b'); %plot contact space

%%
plot(shp_1)
plot(shp_2)
plot(shp_mink)
hold off
alpha(0.5)

%%

prob_max_poly = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2)

exact_prob_x = exact_contact_probability_pure_translation(mu, Sigma, s1, s2, 10000);
%%
figure
 hold on
  s1.PlotShape('b', 0.4);
  s2.PlotShape('g', 0.4);
  %%
X = reshape(x_mink(1,:), N(1), N(2)); %change array to matrix form
Y = reshape(x_mink(2,:), N(1), N(2));
Z = reshape(x_mink(3,:), N(1), N(2));
surf(X, Y, Z,...
 'EdgeColor', 'k', 'EdgeAlpha', 0.2,...
 'FaceAlpha', 0.1, 'FaceColor', 'y'); %plot contact space

%%
n = 100;
mu=s2.tc;
Sigma = zeros(3);
    Sigma(1,1) = 0.0001;
    Sigma(2,2) = 0.0009;
    Sigma(3,3) = 0.0025;
xi_rand = mvnrnd(mu, Sigma, n);
for i = 1:size(xi_rand, 1)
   
    s2.tc = xi_rand(i,1:3)';

    s2.PlotShape('r', 0.1);
end
