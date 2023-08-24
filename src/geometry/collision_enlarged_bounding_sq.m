function flag = collision_enlarged_bounding_sq(s1, s2, mu, Sigma, confidenceLevel)
% This function computes the collision status using enlarged bounding
% convex mesh of s2 based on its pose error propability distribution

bounding_s2 = enlarged_bounding_superquadrics(s2, mu, Sigma, confidenceLevel);
flag=collision_mesh(s1, bounding_s2);

end