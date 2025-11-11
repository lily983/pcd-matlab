function [mf, Sigmaf] = get_bounding_ellip(e1, e2)
% This function gets the bounding ellip for the Minkowski sum of e1 and e2
% e1 and e2 are ellipsoids (3D) or ellipse (2D)

% Check dimension of e1 and e2
dimension = size(e1.a,2);

if dimension==2
    if isequal(e1.eps, 1)==false || isequal(e2.eps, 1)==false
        error("Inputs are not ellipsoid, eps")
    end
    Q1 = angle2rotm(e1.ang) * diag(e1.a).^2 * angle2rotm(e1.ang)';
    Q2 = angle2rotm(e2.ang) * diag(e2.a).^2 * angle2rotm(e2.ang)';
elseif dimension==3
    if isequal(e1.eps, ones(1, 2))==false || isequal(e2.eps, ones(1, 2))==false
        error("Inputs are not ellipsoid, eps")
    end
    Q1 = quat2rotm(e1.q) * diag(e1.a).^2 * quat2rotm(e1.q)';
    Q2 = quat2rotm(e2.q) * diag(e2.a).^2 * quat2rotm(e2.q)';
end

Sigmaf = (1 + sqrt(trace(Q2)/trace(Q1))) * Q1 + (1 + sqrt(trace(Q1)/trace(Q2))) * Q2;
mf = zeros(dimension, 1);
end