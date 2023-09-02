function planningBenchmark
base_path = "/home/xiaoli/pcd_ws/src/cfc_collision_ros/data/result/pcd_planning_benchmark/";
%%
PCDEB95_narrow_large_path = base_path + "PCDEB95_narrow_large.csv";
PCDEB95_narrow_small_path = base_path + "PCDEB95_narrow_small.csv";

PCDEB95_sparse_large_path = base_path + "PCDEB95_sparse_large.csv";
PCDEB95_sparse_small_path = base_path + "PCDEB95_sparse_small.csv";

PCDGMM_narrow_large_path = base_path + "PCDGMM_narrow_large.csv";
PCDGMM_narrow_small_path = base_path + "PCDGMM_narrow_small.csv";

PCDGMM_sparse_large_path = base_path + "PCDGMM_sparse_large.csv";
PCDGMM_sparse_small_path = base_path + "PCDGMM_sparse_small.csv";
%%
PCDEB95_narrow_large = readtable(PCDEB95_narrow_large_path);
PCDEB95_narrow_small = readtable(PCDEB95_narrow_small_path);

PCDEB95_sparse_large = readtable(PCDEB95_sparse_large_path);
PCDEB95_sparse_small = readtable(PCDEB95_sparse_small_path);

PCDGMM_narrow_large = readtable(PCDGMM_narrow_large_path);
PCDGMM_narrow_small = readtable(PCDGMM_narrow_small_path);

PCDGMM_sparse_large = readtable(PCDGMM_sparse_large_path);
PCDGMM_sparse_small = readtable(PCDGMM_sparse_small_path);
%%
PCDEB95_narrow_large_path_length = table2array(PCDEB95_narrow_large(1:size(PCDEB95_narrow_large,1)-2,10));
PCDEB95_narrow_small_path_length = table2array(PCDEB95_narrow_small(1:size(PCDEB95_narrow_small,1)-2,10));

PCDEB95_sparse_large_path_length = table2array(PCDEB95_sparse_large(1:size(PCDEB95_sparse_large,1)-2,10));
PCDEB95_sparse_small_path_length = table2array(PCDEB95_sparse_small(1:size(PCDEB95_sparse_small,1)-2,10));

PCDGMM_narrow_large_path_length = table2array(PCDGMM_narrow_large(1:size(PCDGMM_narrow_large,1)-2,10));
PCDGMM_narrow_small_path_length = table2array(PCDGMM_narrow_small(1:size(PCDGMM_narrow_small,1)-2,10));

PCDGMM_sparse_large_path_length = table2array(PCDGMM_sparse_large(1:size(PCDGMM_sparse_large,1)-2,10));
PCDGMM_sparse_small_path_length = table2array(PCDGMM_sparse_small(1:size(PCDGMM_sparse_small,1)-2,10));
%%
ax1=subplot(2,2,1)
boxplot([PCDEB95_narrow_large_path_length, PCDGMM_narrow_large_path_length],'Labels',{'PCD-EB95','PCD-GMM'});
title('Narrow environment with large inaccurace')

ax2=subplot(2,2,2)
boxplot([PCDEB95_narrow_small_path_length, PCDGMM_narrow_small_path_length],'Labels',{'PCD-EB95','PCD-GMM'});
title('Narrow environment with small inaccurace')


ax3=subplot(2,2,3)
boxplot([PCDEB95_sparse_large_path_length, PCDGMM_sparse_large_path_length],'Labels',{'PCD-EB95','PCD-GMM'});
title('Sparse environment with large inaccurace')


ax4=subplot(2,2,4)
boxplot([PCDEB95_sparse_small_path_length, PCDGMM_sparse_small_path_length],'Labels',{'PCD-EB95','PCD-GMM'});
title('Sparse environment with small inaccurace')

linkaxes([ax1 ax2 ax3 ax4],'xy')

end