gh1=subplot(2,2,1);
hold on;
visualizeResults(resultsSphereSmall, 'sphere');
title("Sphere-Sphere case with small perception uncertainty")
hold off;

gh2=subplot(2,2,2);
hold on;
visualizeResults(resultsEllipSmall, 'ellipsoid');
title("Ellipsoid-Ellipsoid case with small perception uncertainty")
hold off;

gh3=subplot(2,2,3);
hold on;
visualizeResults(resultsSphereLarge, 'sphere');
title("Sphere-Sphere case with large perception uncertainty")
hold off;

gh4=subplot(2,2,4);
hold on;
visualizeResults(resultsEllipLarge, 'ellipsoid');
title("Ellipsoid-Ellipsoid case with large perception uncertainty")
hold off;

linkaxes([gh1 gh2 gh3 gh4], 'xy')
