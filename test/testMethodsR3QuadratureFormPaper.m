clc; clear; close all;
% add_path()

sampleNumber = 10;

%% Generate sampleNumber random samples for spheres
resultsSphere = getResultsR3QuadratureFormPaper(sampleNumber, 'sphere');
visualizeResults(resultsSphere, 'sphere');

%% Generate sampleNumber random samples for ellipsoids
resultsEllip = getResultsR3QuadratureFormPaper(sampleNumber, 'ellipsoid');
visualizeResults(resultsEllip, 'ellipsoid');


 