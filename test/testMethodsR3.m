clc; clear; close all;
% add_path()

sampleNumber = 1e+01;

%% Choose the position error distribution: sparse or concentrate
% distribution = 'small';
distribution = 'large';

%% Generate sampleNumber random samples for spheres
Sigma = positionErrorDistribution(distribution);
resultsSphere = getResultsR3(sampleNumber, 'sphere', Sigma);
visualizeResults(resultsSphere, 'sphere');

%% Generate sampleNumber random samples for ellipsoids
Sigma = positionErrorDistribution(distribution);
resultsEllip = getResultsR3(sampleNumber, 'ellipsoid', Sigma);
visualizeResults(resultsEllip, 'ellipsoid');


