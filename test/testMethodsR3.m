clc; clear; close all;
add_path()

sizeScale = 1;
sampleNumber = 10;

%% Generate 50 random samples for spheres with size range (1, 10)
resultsSphere = getResultsR3(sizeScale, sampleNumber, 'sphere');
visualizeResults(resultsSphere, 'sphere');

%% Generate 50 random samples for ellipsoid with size range (1, 10)
resultsEllip = getResultsR3(sizeScale, sampleNumber, 'ellip');
visualizeResults(resultsEllip, 'ellip');