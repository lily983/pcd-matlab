close all;clear;clc;

N = 20;
sizeScale = 1;

for i=1:N

    quat=rand(1,4);
    quat = quat/norm(quat);

    s1 = SuperQuadrics({(rand(1)+1)*1*sizeScale*ones(1,3), [1,1], [0, 0]...
    rand(3,1)*0.3*sizeScale, quat, [20, 20]});

    quat=rand(1,4);
    quat = quat/norm(quat);

    s2 = SuperQuadrics({(rand(1)+1)*1*sizeScale*ones(1,3),  [1,1], [0, 0]...
       rand(3,1)*3.6*sizeScale, quat,[20,20]});

    [flag, dist, pt_cls, condition] = collision_cfc(s1,s2, 'constrained');
 
    mu=zeros(3,1);
    
    Sigma = zeros(3);
    Sigma(1,1) = 8.0000e-00;
    Sigma(2,2) = 8.0000e-00;
    Sigma(3,3) = 8.0000e-00;
    Sigma = Sigma.*0.01*sizeScale^2;
    
    v =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
    
    try
       [pdf_max_x, x_max] = max_contact_probability_pure_translation(s2.tc, Sigma, s1, s2, flag);
        prob_max_x = v * pdf_max_x;
    catch
        prob_max_x = NaN;
    end
    
    try
        prob_max_poly = max_contact_probability_polyhedron_pure_translation(Sigma, s1, s2);
    catch
        prob_max_poly = NaN;
    end
        
     try
         mu_g = [eye(3), s2.tc; 0, 0, 0, 1];
         Sigma_g = zeros(6,6); Sigma_g(4:6,4:6) = Sigma;
         prob_bounding_sq = collision_enlarged_bounding_sq(s1, s2, mu_g, Sigma_g);
    catch
        prob_bounding_sq = NaN;
     end
     
     try
         [prob_bounding_minksum, weightedNorm] = max_prob_single_gaussian(s1, s2, mu, Sigma);
     catch 
         prob_bounding_minksum = NaN;
     end
     
     try 
         prob_bounding_minksum_GM = max_prob_gaussian_mixture(s1, s2, mu, Sigma, 10, 10, 3);
     catch
         prob_bounding_minksum_GM = NaN;
     end
     
     try 
         prob_ellip_clipped_single_gaussian = max_prob_single_gaussian_clipped_integration(s1, s2, mu, Sigma);
     catch
         prob_ellip_clipped_single_gaussian = NaN;
     end
     
     try 
         prob_ellip_clipped_GM = max_prob_gaussian_mixture_clipped_integration(s1, s2, mu, Sigma);
     catch
         prob_ellip_clipped_GM = NaN;
     end
    
    try
        exact_prob_x = exact_prob_translation(s1, s2, Sigma, 1e+03);
    catch
        exact_prob_x = NaN;
    end
        
    j=1; result(j, i) = round(flag);
    j=j+1; result(j, i) = dist;
    
    j=j+1; result(j, i) = exact_prob_x;
    j=j+1; result(j, i) = prob_max_x;
    j=j+1; result(j, i) = prob_max_poly;
    j=j+1; result(j, i) = prob_bounding_sq; 
    j=j+1; result(j, i) = prob_bounding_minksum;
    j=j+1; result(j, i) = prob_bounding_minksum_GM;
    j=j+1; result(j, i) = weightedNorm;

    j=j+2; result(j, i) = s1.a(1);
    j=j+1; result(j:j+2, i) = s1.tc;
    j=j+3; result(j:j+3, i) =s1.q;
    j=j+5; result(j, i) = s2.a(1);
    j=j+1; result(j:j+2, i) = s2.tc;
    j=j+3; result(j:j+3, i) =s2.q;
    
    j=j+5; result(j, i) = prob_ellip_clipped_single_gaussian;
    j=j+1; result(j, i) = prob_ellip_clipped_GM;

end   

%% Plot results
%  N=500;
flag = result(1, 1:N);
dist = result(2, 1:N);
exact = result(3, 1:N);
x_max = result(4, 1:N);
convex_hull = result(5, 1:N);
enlarged_bounding_sq = result(6,1:N);
our_bound = result(7, 1:N);
our_bound_GMM = result(8, 1:N);
weighted_norm = result(9, 1:N);
our_bound_EC = result(29, 1:N);
our_bound_GMM_EC=result(30, 1:N);

% 
% for i=1:N
%     if (convex_hull(i)>1)
%         convex_hull(i)=1;
%     end
%     if(x_max(i)>1)
%         x_max(i)=1;
%     end
%     if(our_bound_GMM(i)>1)
%         our_bound_GMM(i)=1;
%     end
%     if(our_bound(i)>1)
%         our_bound(i)=1;
%     end
%     if(our_bound_GMM_EC(i)>1)
%         our_bound_GMM_EC(i)=1;
%     end
%     if(our_bound_EC(i)>1)
%         our_bound_EC(i)=1;
%     end
% end

% %% sort by weighted norm
% [sort_norm, index] = sort(weighted_norm, 2, 'descend');
[sort_exact, index] = sort(exact, 2);

for i=1:N
    sort_flag(i) = flag(index(i));
    sort_exact(i) = exact(index(i));
    sort_our_bound(i) = our_bound(index(i));
    sort_our_bound_GMM(i) = our_bound_GMM(index(i));
    sort_x_max(i) = x_max(index(i));
    sort_convex_hull(i) = convex_hull(index(i));
    sort_enlarged_sq(i) = enlarged_bounding_sq(index(i));
    sort_our_bound_EC(i) = our_bound_EC(index(i));
    sort_our_bound_GMM_EC(i) = our_bound_GMM_EC(index(i));
end

%%
figure; hold on;
plot(1:N, sort_exact(1:N), 'k')
% plot(1:N, sort_convex_hull(1:N), 'g')
% plot(1:N, sort_enlarged_sq(1:N), 'y')
% plot(1:N,  sort_x_max(1:N), 'm')
plot(1:N, sort_our_bound(1:N), 'r')
plot(1:N, sort_our_bound_GMM(1:N), 'b')
% plot(1:N,  sort_norm(1:N), 'm');
% plot(1:N, sort_flag(1:N),'c')
% plot(1:N, ones(1, N),'y')
plot(1:N, sort_our_bound_EC(1:N), 'g')
plot(1:N, sort_our_bound_GMM_EC(1:N), 'y')

%% 
figure; hold on;
scatter(1:N, sort_exact(1:N),...
     'MarkerFaceColor', hex2rgb('073b4c')./255,...
     'AlphaData', sort_exact(1:N),...
     'MarkerFaceAlpha', 'flat',...
     'MarkerEdgeAlpha', 0.0);

scatter(1:N, sort_enlarged_sq(1:N),...
     'MarkerFaceColor', hex2rgb('ffd166')./255,...
     'AlphaData', sort_enlarged_sq(1:N),...
     'MarkerFaceAlpha', 'flat',...
     'MarkerEdgeAlpha', 0.0);
 
scatter(1:N, sort_convex_hull(1:N),...
     'MarkerFaceColor', hex2rgb('06d6a0')./255,...
     'AlphaData', sort_convex_hull(1:N),...
     'MarkerFaceAlpha', 'flat',...
     'MarkerEdgeAlpha', 0.0);
 
scatter(1:N, sort_our_bound(1:N),...
     'MarkerFaceColor', hex2rgb('ef476f')./255,...
     'AlphaData', sort_our_bound(1:N),...
     'MarkerFaceAlpha', 'flat',...
     'MarkerEdgeAlpha', 0.0);
 
scatter(1:N, sort_our_bound_GMM(1:N),...
     'MarkerFaceColor', hex2rgb('118ab2 ')./255,...
     'AlphaData', sort_our_bound_GMM(1:N),...
     'MarkerFaceAlpha', 'flat',...
     'MarkerEdgeAlpha', 0.0);
%  %%
%  figure; hold on;
 markerAlpha = 0.4;
 markerSize = 30;
scatter(1:N, sort_exact(1:N),...
    markerSize,...
     'MarkerFaceColor', hex2rgb('073b4c')./255,...
     'MarkerFaceAlpha',markerAlpha,...
     'MarkerEdgeAlpha', 0.0);

scatter(1:N, sort_enlarged_sq(1:N),...
     markerSize,...
     'MarkerFaceColor', hex2rgb('fb8500')./255,...
     'MarkerFaceAlpha', markerAlpha,...
     'MarkerEdgeAlpha', 0.0);
 
scatter(1:N, sort_convex_hull(1:N),...
    markerSize,...
     'MarkerFaceColor', hex2rgb('06d6a0')./255,...
     'MarkerFaceAlpha',markerAlpha,...
     'MarkerEdgeAlpha', 0.0);
 
scatter(1:N, sort_our_bound(1:N),...
    markerSize,...
     'MarkerFaceColor', hex2rgb('ef476f')./255,...
     'MarkerFaceAlpha', markerAlpha,...
     'MarkerEdgeAlpha', 0.0);
 
scatter(1:N, sort_our_bound_GMM(1:N),...
    markerSize,...
     'MarkerFaceColor', hex2rgb('118ab2 ')./255,...
     'MarkerFaceAlpha', markerAlpha,...
     'MarkerEdgeAlpha', 0.0);


% figure; hold on;
VisualizeFitResult(sort_exact, hex2rgb('073b4c')./255)
% VisualizeFitResult(sort_enlarged_sq,hex2rgb('ffd166')./255)
% plot(1:N, sort_enlarged_sq(1:N), 'Color', hex2rgb('ffd166')./255, 'LineWidth', 2)
VisualizeFitResult(sort_convex_hull, hex2rgb('06d6a0')./255)
% VisualizeFitResult(sort_x_max, 'm')
VisualizeFitResult(sort_our_bound, hex2rgb('ef476f')./255)
VisualizeFitResult(sort_our_bound_GMM, hex2rgb('118ab2 ')./255)
% VisualizeFitResult(sort_enlarged_sq, 'y')

%%
% figure; hold on;
% VisualizeFitResult(sort_convex_hull, hex2rgb('06d6a0')./255)
% Define the "shading"
% Note how each x_points(i) corresponds to y_points(i)

Sigma = zeros(3);
Sigma(1,1) = 8.0000e-00;
Sigma(2,2) = 8.0000e-00;
Sigma(3,3) = 8.0000e-00;
Sigma = Sigma.*0.01*sizeScale^2;

%size 0.1: n=137;n2=238;
%size 1: n=151; n2=262;
n=137;
prob=sort_exact(n);
n2=238;
prob2=sort_exact(n2);

x1_points = [0, 0, n, n];  
y1_points = [0, prob, prob, 0];
color1 = hex2rgb('06d6a0')./255;
a1 = fill(x1_points, y1_points, color1);
a1.FaceAlpha = 0.6;
a1.EdgeAlpha = 0.0;

x2_points = [n, n, n2, n2];  
y2_points = [prob, 1, 1, prob];
color2 =  hex2rgb('118ab2 ')./255;
a2 = fill(x2_points, y2_points, color2);
a2.FaceAlpha = 0.2;
a2.EdgeAlpha = 0.0;
% 
vertical_line = linspace(0, 1, 100);
plot(n*ones(size(vertical_line)), vertical_line, 'Color', hex2rgb('073b4c')./255);

horizontal_line = prob * ones(1, n2);
plot(1:n2, horizontal_line, 'Color', hex2rgb('073b4c')./255);

vertical_line2 = linspace(prob, 1, 100);
plot(n2*ones(size(vertical_line)), vertical_line2, 'Color', hex2rgb('073b4c')./255);

horizontal_line2 = prob2 * ones(1, n2);
plot(1:n2, horizontal_line2, 'Color', hex2rgb('073b4c')./255, 'LineStyle',":", 'LineWidth', 2);

horizontal_line3 = prob * ones(1, N-n2+1);
plot(n2:N, horizontal_line3, 'Color', hex2rgb('073b4c')./255, 'LineStyle',":", 'LineWidth', 2);

horizontal_line4 = prob2 * ones(1, N-n2+1);
plot(n2:N, horizontal_line4, 'Color', hex2rgb('073b4c')./255, 'LineStyle',":", 'LineWidth', 2);

yticks([0 prob  0.2  0.4 prob2  0.6  0.8  1.0])
ylim([0 1])
xlim([1 500])
