clc; clear; close all;
add_path()

sampleNumber = 1e+02;

%%
% % Mute key-value pair if don't want to test any method
% % key is the name of PCD method
% % value is the color of PCD method
M = containers.Map('KeyType','char','ValueType','any');
M('Exact') = hex2rgb('000000'); % dark blue
M('Exact_two') = hex2rgb('000000'); % black 
M('Fast_sampling') = hex2rgb('f59494') ;% pink
M('Divergence_mesh') = hex2rgb('45AC59'); % light green
M('Maxpdf') = hex2rgb('f58612'); % origin
M('EB_99') = hex2rgb('9f580f'); % brown
M('EB_95') = hex2rgb('98646b'); % grape
M('Quadratic_bound') = hex2rgb('8f90dd') ;% light purple

% % Methods: LCC_center_point, LCC_tangent_point, LCC_closed_point, Maxpdf_SQ
% % are extended methods from previous works by us
M('LCC_center_point') = hex2rgb('2023c7'); % bright purple
M('LCC_tangent_point') = hex2rgb('8a8686'); % light grey
M('LCC_closed_point') = hex2rgb('EBBF00'); % ginger yellow
M('LCC_center_point_cfc') = hex2rgb('1d5406'); % dark green
% M('Maxpdf_SQ') = hex2rgb('4996db'); % light blue

% % GMF is the closed-form solution proposed by us
M('GMF') = hex2rgb('ff0010') ;% red

%% Choose the position error distribution: 'none', 'small', 'large'
distribution = 'small';
% distribution = 'large';
% distribution = 'none';
Sigma1 = positionErrorDistribution('small');
Sigma2 = positionErrorDistribution('large');

TwoErrorCase = 0;
if all(diag(Sigma1))
    TwoErrorCase = 1;
end

%% Generate sampleNumber random samples for spheres
resultsSphere = getResultsR3(sampleNumber, 'sphere', Sigma1, Sigma2, M.keys);

visualizeResults(resultsSphere, M.keys, M.values, sampleNumber);

[ratio_sphere, time_sphere, ub_sphere] = calculatePCDBenchmark(resultsSphere, M.keys, TwoErrorCase)

%% Generate sampleNumber random samples for ellipsoids
resultsEllip = getResultsR3(sampleNumber, 'ellipsoid', Sigma1, Sigma2, M.keys);

visualizeResults(resultsEllip, M.keys, M.values, sampleNumber);

[ratio_ellip, time_ellip, ub_ellip] = calculatePCDBenchmark(resultsEllip, M.keys, TwoErrorCase)

%% Generate sampleNumber random samples for SQ
resultsSQ = getResultsR3(sampleNumber, 'superquadrics', Sigma1, Sigma2, M.keys);

visualizeResults(resultsSQ, M.keys, M.values, sampleNumber);

[ratio_SQ, time_SQ, ub_SQ] = calculatePCDBenchmark(resultsSQ, M.keys, TwoErrorCase)