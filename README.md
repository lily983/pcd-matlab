# prob-collision
MATLAB implementation of the proposed probabilistic collision detection (PCD) algorithm using Principal Kinematic Formula and Kinematic Inequality.

## Dependencies
- Manifold optimization toolbox: [manopt 7.0.0](https://www.manopt.org/index.html)

## Test scripts
- [Bug] [test_prob_collision.m](test/test_prob_collision.m): Test for the proposed PCD algorithm
- [test_poly3d_integral_mean_curvature.m](test/test_poly3d_integral_mean_curvature.m): Test for integral of mean curvature for 3D polyhedra
- [main_pkf_3d.m](test/main_pkf_3d.m): Main script for computing Principal Kinematic Formula (PKF) in 3D