function [pkf, s1_geom, s2_geom] = pkf_3d(s1, s2, ispoly, use_fcl_ellipsoid)

%% Volume, surface area and integral of mean curvature
if ~ispoly
    s1.SetNumVertices([200,200]);
    s2.SetNumVertices([200,200]);
end

% Generate polyhedra for superquadric and ellipsoid objects
s1_poly = s1.GetPoints()';
if use_fcl_ellipsoid
    s2_poly = get_fcl_ellipsoid(s2);
else
    s2_poly = s2.GetPoints()';
end

% Volume
if ~ispoly  
    s1_geom.volume = s1.GetVolume();
    s2_geom.volume = s2.GetVolume();
    k1 = convhull(s1_poly);
    k2 = convhull(s2_poly);
else
    [k1, s1_geom.volume] = convhull(s1_poly);
    [k2, s2_geom.volume] = convhull(s2_poly);
end

% Surface area
s1_geom.area = poly3d_surface_area(s1_poly, k1);
s2_geom.area = poly3d_surface_area(s2_poly, k2);

% Integral of mean curvature
s1_geom.int_curvature = poly3d_integral_mean_curvature(s1_poly, k1);
s2_geom.int_curvature = poly3d_integral_mean_curvature(s2_poly, k2);

%% Compute PKF
pkf = 8 * pi^2 * (s1_geom.volume + s2_geom.volume) +...
    2 * pi * (s1_geom.area * s2_geom.int_curvature +...
    s2_geom.area * s1_geom.int_curvature);