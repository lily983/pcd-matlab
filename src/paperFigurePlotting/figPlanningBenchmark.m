function planningBenchmark
base_path = "/home/xiaoli/pcd_ws/src/cfc_collision_ros/data/result/pcd_planning_benchmark/";
%%
PCDEB95_narrow_large_path = base_path + "PCDEB95_narrow_large.csv";
PCDEB95_narrow_small_path = base_path + "PCDEB95_narrow_small.csv";
% 
% PCDEB95_sparse_large_path = base_path + "PCDEB95_sparse_large.csv";
% PCDEB95_sparse_small_path = base_path + "PCDEB95_sparse_small.csv";

PCDEB95_clamp_large_path = base_path + "PCDEB95_clamp_large.csv";
PCDEB95_clamp_small_path = base_path + "PCDEB95_clamp_small.csv";

PCDGMM_narrow_large_path = base_path + "PCDGMM_narrow_large.csv";
PCDGMM_narrow_small_path = base_path + "PCDGMM_narrow_small.csv";

% PCDGMM_sparse_large_path = base_path + "PCDGMM_sparse_large.csv";
% PCDGMM_sparse_small_path = base_path + "PCDGMM_sparse_small.csv";

PCDGMM_clamp_large_path = base_path + "PCDGMM_clamp_large.csv";
PCDGMM_clamp_small_path = base_path + "PCDGMM_clamp_small.csv";
%%
PCDEB95_narrow_large = readtable(PCDEB95_narrow_large_path);
PCDEB95_narrow_small = readtable(PCDEB95_narrow_small_path);

% PCDEB95_sparse_large = readtable(PCDEB95_sparse_large_path);
% PCDEB95_sparse_small = readtable(PCDEB95_sparse_small_path);

PCDEB95_clamp_large = readtable(PCDEB95_clamp_large_path);
PCDEB95_clamp_small = readtable(PCDEB95_clamp_small_path);

PCDGMM_narrow_large = readtable(PCDGMM_narrow_large_path);
PCDGMM_narrow_small = readtable(PCDGMM_narrow_small_path);

% PCDGMM_sparse_large = readtable(PCDGMM_sparse_large_path);
% PCDGMM_sparse_small = readtable(PCDGMM_sparse_small_path);

PCDGMM_clamp_large = readtable(PCDGMM_clamp_large_path);
PCDGMM_clamp_small = readtable(PCDGMM_clamp_small_path);
%%
PCDEB95_narrow_large_path_length = table2array(PCDEB95_narrow_large(1:size(PCDEB95_narrow_large,1)-4,10));
PCDEB95_narrow_small_path_length = table2array(PCDEB95_narrow_small(1:size(PCDEB95_narrow_small,1)-4,10));

% PCDEB95_sparse_large_path_length = table2array(PCDEB95_sparse_large(1:size(PCDEB95_sparse_large,1)-4,10));
% PCDEB95_sparse_small_path_length = table2array(PCDEB95_sparse_small(1:size(PCDEB95_sparse_small,1)-4,10));

PCDEB95_clamp_large_path_length = table2array(PCDEB95_clamp_large(1:size(PCDEB95_clamp_large,1)-4,10));
PCDEB95_clamp_small_path_length = table2array(PCDEB95_clamp_small(1:size(PCDEB95_clamp_small,1)-4,10));

PCDGMM_narrow_large_path_length = table2array(PCDGMM_narrow_large(1:size(PCDGMM_narrow_large,1)-4,10));
PCDGMM_narrow_small_path_length = table2array(PCDGMM_narrow_small(1:size(PCDGMM_narrow_small,1)-4,10));

% PCDGMM_sparse_large_path_length = table2array(PCDGMM_sparse_large(1:size(PCDGMM_sparse_large,1)-4,10));
% PCDGMM_sparse_small_path_length = table2array(PCDGMM_sparse_small(1:size(PCDGMM_sparse_small,1)-4,10));

PCDGMM_clamp_large_path_length = table2array(PCDGMM_clamp_large(1:size(PCDGMM_clamp_large,1)-4,10));
PCDGMM_clamp_small_path_length = table2array(PCDGMM_clamp_small(1:size(PCDGMM_clamp_small,1)-4,10));
%%
figure;
fontSize=8;
colors = validatecolor(["#0064B0" "#d14f53"], 'multiple');

ax1=subplot(2,2,1)
boxplot([PCDEB95_narrow_large_path_length, PCDGMM_narrow_large_path_length],'Labels',{'PCD-EB95','PCD-GMM'},...
    'BoxStyle','outline','Colors',colors);
set(ax1, 'FontName', 'Times', 'FontSize', fontSize);
ylabel('Average path length')
% title('Narrow environment with large environment perception uncertainty');
title('Narrow environment (large)')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
end

ax2=subplot(2,2,2)
boxplot([PCDEB95_narrow_small_path_length, PCDGMM_narrow_small_path_length],'Labels',{'PCD-EB95','PCD-GMM'},...
    'BoxStyle','outline','Colors',colors);
% plotTable(PCDEB95_narrow_small_path_length,PCDGMM_narrow_small_path_length);
set(ax2, 'FontName', 'Times', 'FontSize', fontSize);
ylabel('Average path length')
% title('Narrow environment with small environment perception uncertainty')
title('Narrow environment (small)')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
end

% ax3=subplot(3,2,3)
% boxplot([PCDEB95_sparse_large_path_length, PCDGMM_sparse_large_path_length],'Labels',{'PCD-EB95','PCD-GMM'},...
%     'BoxStyle','outline','Colors',colors);
% % plotTable(PCDEB95_sparse_large_path_length,PCDGMM_sparse_large_path_length);
% set(ax3, 'FontName', 'Times', 'FontSize', fontSize);
% ylabel('Average path length')
% title('Sparse environment with large environment perception uncertainty')
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
% end


% ax4=subplot(3,2,4)
% boxplot([PCDEB95_sparse_small_path_length, PCDGMM_sparse_small_path_length],'Labels',{'PCD-EB95','PCD-GMM'},...
%     'BoxStyle','outline','Colors',colors);
% % plotTable(PCDEB95_sparse_small_path_length,PCDGMM_sparse_small_path_length);
% set(ax4, 'FontName', 'Times', 'FontSize', fontSize);
% ylabel('Average path length')
% title('Sparse environment with small perception uncertainty')
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
% end

ax5=subplot(2,2,3)
boxplot([PCDEB95_clamp_large_path_length, PCDGMM_clamp_large_path_length],'Labels',{'PCD-EB95','PCD-GMM'},...
    'BoxStyle','outline','Colors',colors);
% plotTable(PCDEB95_clamp_large_path_length,PCDGMM_clamp_large_path_length);
set(ax5, 'FontName', 'Times', 'FontSize', fontSize);
ylabel('Average path length')
% title('Clamp environment with large environment perception uncertainty')
title('Clamped environment (large)')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
end

ax6=subplot(2,2,4)
boxplot([PCDEB95_clamp_small_path_length, PCDGMM_clamp_small_path_length],'Labels',{'PCD-EB95','PCD-GMM'},...
    'BoxStyle','outline','Colors',colors);
% plotTable(PCDEB95_clamp_small_path_length,PCDGMM_clamp_small_path_length);
set(ax6, 'FontName', 'Times', 'FontSize', fontSize);
ylabel('Average path length')
% title('Clamp environment with small environment perception uncertainty')
title('Clamped environment (small)')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
end

linkaxes([ax1 ax2  ax5 ax6],'xy')

end

function box_plot=plotTable(PCDEB95, PCDGMM)
plot_table = table;

plot_table.Name=[repmat("PCD-EB95", 100, 1); repmat("PCD-GMM", 100, 1)];
plot_table.Data=[PCDEB95; PCDGMM];

name=categorical(plot_table.Name);
box_plot = boxchart(name, plot_table.Data);
set(get(get(box_plot, 'Parent'), 'yLabel'), 'String', 'Path length',  'FontName', 'Times New Roman', 'FontSize', 10);
% box_plot(2).SeriesIndex = 2;
end