clc; clear; close all;
% add_path()

sampleNumber = 1e+1;

keySet = {'Exact', 'LCC_center_point', 'LCC_closed_point', 'LCC_tangent_point'};

colorSet = {hex2rgb('45498C'),  hex2rgb('f59494'), hex2rgb('45AC59'),  hex2rgb('EBBF00')};

%% Choose the position error distribution: sparse or concentrate
% distribution = 'small';
distribution = 'large';

Sigma = positionErrorDistribution(distribution);
%% Generate sampleNumber random samples for spheres
% Sigma = positionErrorDistribution(distribution);
resultsSQ = getResultsR3(sampleNumber, 'superquadrics', Sigma, keySet);

visualizeResults(resultsSQ, 'superquadrics', keySet, colorSet, sampleNumber);

% %% Generate sampleNumber random samples for spheres
% % Sigma = positionErrorDistribution(distribution);
% resultsSphere = getResultsR3(sampleNumber, 'sphere', Sigma);
% visualizeResults(resultsSphere, 'sphere');
% 
% %% Generate sampleNumber random samples for ellipsoids
% % Sigma = positionErrorDistribution(distribution);
% resultsEllip = getResultsR3(sampleNumber, 'ellipsoid', Sigma);
% visualizeResults(resultsEllip, 'ellipsoid');


