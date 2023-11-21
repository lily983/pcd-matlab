function [mf, Sigmaf] = get_bounding_ellip(s1, s2)
% This function gets the bounding ellip for the Minkowski sum of s1 and s2
% s1 and s2 are ellipsoids (3D) or ellipse (2D)

% Check dimension of s1 and s2
dimension = size(s1.a,2);

if dimension==2
    if isequal(s1.eps, 1)==false || isequal(s2.eps, 1)==false
        error("Inputs are not ellipsoid, eps")
    end
    Q1 = angle2rotm(s1.ang) * diag(s1.a).^2 * angle2rotm(s1.ang)';
    Q2 = angle2rotm(s2.ang) * diag(s2.a).^2 * angle2rotm(s2.ang)';
elseif dimension==3
    if isequal(s1.eps, ones(1, 2))==false || isequal(s2.eps, ones(1, 2))==false
        error("Inputs are not ellipsoid, eps")
    end
    Q1 = quat2rotm(s1.q) * diag(s1.a).^2 * quat2rotm(s1.q)';
    Q2 = quat2rotm(s2.q) * diag(s2.a).^2 * quat2rotm(s2.q)';
end

Sigmaf = (1 + sqrt(trace(Q2)/trace(Q1))) * Q1 + (1 + sqrt(trace(Q1)/trace(Q2))) * Q2;
mf = zeros(dimension, 1);
end