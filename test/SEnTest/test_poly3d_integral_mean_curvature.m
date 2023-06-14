close all; clear; clc;
add_path()

% Test file for the integral of mean curvature for convex polyhedron

% N^2: number of vertices on the surface
N = [100,100];

%% Case 1: Sphere
r = 10*rand;
s1 = SuperQuadrics({r*ones(1,3), [1,1], [0,0], zeros(3,1), [1,0,0,0], N});

pts1 = s1.GetPoints()';
k1 = convhull(pts1);

int_curv1 = poly3d_integral_mean_curvature(pts1, k1)
int_curv1_true = 4*pi*r

%% Case 2: Convex cylinder
r2 = 10*rand;
h2_half = 20*rand;
s2 = SuperQuadrics({[r2,r2,h2_half], [.1,1], [0,0], zeros(3,1), [1,0,0,0], N});

pts2 = s2.GetPoints()';
k2 = convhull(pts2);

int_curv2 = poly3d_integral_mean_curvature(pts2, k2)
int_curv2_true = pi*(pi*r2+2*h2_half)

%% Case 3: Ellipsoid of revolution
% 0 < gamma < 1 
a3 = 10*rand;
gamma3 = rand;
s3 = SuperQuadrics({[a3,a3,gamma3*a3], [1,1], [0,0], zeros(3,1), [1,0,0,0], N});

pts3 = s3.GetPoints()';
k3 = convhull(pts3);

int_curv3 = poly3d_integral_mean_curvature(pts3, k3)
int_curv3_true = 2*pi*(gamma3+(1-gamma3^2)^(-1/2)*acos(gamma3))*a3

% gamma > 1
a4 = 10*rand;
gamma4 = 20*rand+1;
s4 = SuperQuadrics({[a4,a4,gamma4*a4], [1,1], [0,0], zeros(3,1), [1,0,0,0], N});

pts4 = s4.GetPoints()';
k4 = convhull(pts4);

int_curv4 = poly3d_integral_mean_curvature(pts4, k4)
int_curv4_true = 2*pi*(gamma4+...
    (gamma4^2-1)^(-1/2)*log(gamma4+(gamma4^2-1)^(1/2)))*a4