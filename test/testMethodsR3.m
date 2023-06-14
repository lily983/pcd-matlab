clc; clear; close all;
add_path()
%% Generate 500 random samples for spheres with size range (1, 10)
resultsSphere = getResultsR3(1, 10, 'sphere');
visualizeResults(resultsSphere, 'sphere');

%% Generate 500 random samples for ellipsoid with size range (1, 10)
resultsEllip = getResultsR3(1, 10, 'ellip');
visualizeResults(resultsEllip, 'ellip');