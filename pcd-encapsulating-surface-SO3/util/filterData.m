function [g_filter, filtered_indices] = filterData(g, mu_R)
% IQR filtering: it doesn't assume the data is Gaussian distributed

R = g(1:3,1:3,:);

d_R = zeros(1, size(R,3));
for i=1:size(R,3)
    d_R(i) = norm(skew2vec(logm_SO(R(:,:,i) \ mu_R)));
end

% Step 1: Calculate the first quartile (Q1), third quartile (Q3), and the interquartile range (IQR)
Q1 = prctile(d_R, 25);  % 25th percentile (first quartile)
Q3 = prctile(d_R, 75);  % 75th percentile (third quartile)
IQR_value = Q3 - Q1;     % Interquartile range

% Step 2: Calculate whisker boundaries (lower and upper bounds)
lower_bound = Q1 - 1.5 * IQR_value;
upper_bound = Q3 + 1.5 * IQR_value;

% Step 3: Filter out data points that are outside the whisker boundaries
% outlier_indices = find(d_g < lower_bound | d_g > upper_bound); % Indices of outliers
% outliers = d_g(outlier_indices); % The actual outlier values

filtered_indices = d_R >= lower_bound & d_R <= upper_bound; 
% filtered = d_g(filtered_indices);

g_filter = g(:,:,filtered_indices);

end