% Test script for fitting approximation errors and deficit ratios using
% methods 1-10
%  Author:
%    Xiaoli Wang, wangxiaoli@u.nus.edu

close all;  clc;
add_path();

%% Loading rotations from samples
% T = csv2group(marker38new, 'PCG3'); %T_low
% T = csv2group(summaryposedata, 'PCG3'); %T_medium
T = csv2group(tinychair, 'PCG3'); %T_high

T_filter = filterData(T, eye(3));

R = T_filter(1:3,1:3,:);
mu_SO3 = get_mean_cov(R, 'SO', 1);

%% Loading obj parameters
% axes=[0.2  0.020  0.12707044]' + rand(3,1)*0.1;
axes=rand(3,1) + 0.01*ones(3,1);
axes = axes/2;

ellip = SuperQuadrics({axes, [1,1], [0,0], zeros(3,1), rotm2quat(mu_SO3), [20,20]});

%% Visualize the obj 
figure; hold on; axis equal

for i=1:size(R,3)
    ellip_i = SuperQuadrics({ellip.a, ellip.eps, [0,0], zeros(3,1), rotm2quat(R(:,:,i)), [20,20]});
    ellip_i.PlotShape('b', 0.05,0.3);
end
ellip.PlotShape('y', 0.4,0.2);
%% Generate surface points of ellip at errorneous orientations
ellip_points_cell = {};

for i=1:size(R,3)
    ellip_i = SuperQuadrics({ellip.a, ellip.eps, [0,0], zeros(3,1), rotm2quat(R(:,:,i)), [20,20]});
    ellip_points_cell{i} = ellip_i.GetPoints();
end

%% Range of k values
k_values = 1.0:0.05:1.4;   % from 1.0 to 1.4 with step 0.05
num_k = numel(k_values);

num_methods = 10; 

% Preallocate storage: methods × k
fitting_scores = zeros(num_methods, num_k);
deficit_scores = zeros(num_methods, num_k);

% Loop over k
for ki = 1:num_k
    k = k_values(ki);

    % --- Generate surfaces for all methods at this k ---
    x_new_method = cell(1, num_methods);
    x_new_method{1}  = method1(ellip, R, k);
    x_new_method{2}  = method2(ellip, R, k);
    x_new_method{3}  = method3(ellip, R, k);
    x_new_method{4}  = method4(ellip, R, k);
    x_new_method{5}  = method5(ellip, R, k);
    x_new_method{6}  = method6(ellip, k);
    x_new_method{7}  = method7(ellip, R);
    x_new_method{8}  = method8(ellip, R);
    x_new_method{9}  = method9(ellip, R, k);
    x_new_method{10} = method10(ellip, R, k);

    % --- Compute scores ---
    for m = 1:num_methods
        fitting_scores(m, ki) = metric_fitting_multi_Strue(x_new_method{m}, ellip_points_cell);
        deficit_scores(m, ki) = metric_deficit_multi_Strue(x_new_method{m}, ellip_points_cell);
    end
end

%% ---- Labels + quick sanity ----
method_names = "Method " + string(1:num_methods);
k_labels = compose('k=%.2f', k_values);

% Clamp numerical noise if any metric has tiny negatives / >1 due to voxelization
fitting_scores = max(0, min(1, fitting_scores));
deficit_scores = max(-Inf, min(1, deficit_scores));  % r_def can be <1 but should not exceed 1
%% ---- Overall summary (mean ± std over k) ----
iou_mean = mean(fitting_scores, 2, 'omitnan');
iou_std  = std( fitting_scores, 0, 2, 'omitnan');

rdef_mean = mean(deficit_scores, 2, 'omitnan');
rdef_std  = std( deficit_scores, 0, 2, 'omitnan');

T_overall = table(method_names', iou_mean, iou_std, rdef_mean, rdef_std, ...
    'VariableNames', {'Method','IoU_mean','IoU_std','rdef_mean','rdef_std'});
disp('=== Overall Summary (mean ± std across k) ===');
disp(T_overall);

% Save CSV for thesis tables
writetable(T_overall, 'overall_summary_k_sweep.csv');
%% ---- Heatmaps: methods × k ----
figure('Name','IoU Heatmap (methods × k)','Color','w');
imagesc(fitting_scores, [0 1]); axis tight; colorbar; colormap(parula(numel(colormap)));   % parula is very close to viridis
set(gca,'XTick',1:numel(k_labels),'XTickLabel',k_labels,'XTickLabelRotation',45);
set(gca,'YTick',1:num_methods,'YTickLabel',method_names);
xlabel('k'); ylabel('Method'); title('IoU');
set(gcf,'Position',[100 100 1200 420]);
exportgraphics(gcf,'heatmap_iou_methods_by_k.png','Resolution',300);

figure('Name','r_{def} Heatmap (methods × k)','Color','w');
imagesc(deficit_scores, [min(deficit_scores(:), [], 'omitnan') 1]); axis tight; colorbar;colormap(parula(numel(colormap))); 
set(gca,'XTick',1:numel(k_labels),'XTickLabel',k_labels,'XTickLabelRotation',45);
set(gca,'YTick',1:num_methods,'YTickLabel',method_names);
xlabel('k'); ylabel('Method'); title('Deficit ratio r_{def}');
set(gcf,'Position',[100 100 1200 420]);
exportgraphics(gcf,'heatmap_rdef_methods_by_k.png','Resolution',300);

% CSV exports for heatmap tables
T_iou_grid = array2table(fitting_scores, 'VariableNames', cellstr(k_labels), 'RowNames', cellstr(method_names));
T_rdef_grid = array2table(deficit_scores, 'VariableNames', cellstr(k_labels), 'RowNames', cellstr(method_names));
writetable(T_iou_grid, 'iou_grid_methods_by_k.csv','WriteRowNames',true);
writetable(T_rdef_grid,'rdef_grid_methods_by_k.csv','WriteRowNames',true);

