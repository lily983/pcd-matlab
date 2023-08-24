clc; clear; close all;
% add_path()

sampleNumber = 50;

% Choose the position error distribution: sparse or concentrate
distribution = 'small';
% distribution = 'large';

%% Generate sampleNumber random samples for spheres
Sigma = positionErrorDistribution('sphere', distribution);
resultsSphere = getResultsR3(sampleNumber, 'sphere', Sigma);
visualizeResults(resultsSphere, 'sphere');

%% Generate sampleNumber random samples for ellipsoids
Sigma = positionErrorDistribution('ellipsoid', distribution);
resultsEllip = getResultsR3(sampleNumber, 'ellipsoid', Sigma);
visualizeResults(resultsEllip, 'ellipsoid');


