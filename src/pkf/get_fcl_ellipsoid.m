function vtx = get_fcl_ellipsoid(ellipsoid)
%GET_FCL_ELLIPSOID Computes bounding volume for ellipsoid that is
%implemented in FCL
%
%  Input:
%    ellipsoid: a SuperQuadrics class with eps = [1,1]
%
%  Output:
%    vtx      : 12x3 array, each row represents a vertex on the bounding
%               volume
%
%  Rewritten from FCL library (source: https://flexible-collision-library.github.io/d0/d9e/ellipsoid-inl_8h_source.html)

if ellipsoid.eps(1) ~= 1 || ellipsoid.eps(2) ~= 1
    error('Not an ellipsoid')
end

vtx = nan(12,3);

phi = (1.0 + sqrt(5.0)) / 2;
a = sqrt(3) / (phi * phi);
b = phi * a;

A = ellipsoid.a(1);
B = ellipsoid.a(2);
C = ellipsoid.a(3);

Aa = A * a;
Ab = A * b;
Ba = B * a;
Bb = B * b;
Ca = C * a;
Cb = C * b;

vtx(1,:) = [0, Ba, Cb];
vtx(2,:) = [0, -Ba, Cb];
vtx(3,:) = [0, Ba, -Cb];
vtx(4,:) = [0, -Ba, -Cb];
vtx(5,:) = [Aa, Bb, 0];
vtx(6,:) = [-Aa, Bb, 0];
vtx(7,:) = [Aa, -Bb, 0];
vtx(8,:) = [-Aa, -Bb, 0];
vtx(9,:) = [Ab, 0, Ca];
vtx(10,:) = [Ab, 0, -Ca];
vtx(11,:) = [-Ab, 0, Ca];
vtx(12,:) = [-Ab, 0, -Ca];

vtx = quat2rotm(ellipsoid.q) * vtx' + ellipsoid.tc;
vtx = vtx';
end

