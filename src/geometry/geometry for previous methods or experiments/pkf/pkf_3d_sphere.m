function [pkf, s1_geom, s2_geom] = pkf_3d_sphere(s1, s2)
%{
pkf_3d_sphere: Calculate PKF value for two sphere objects

Inputs:
s1, s2: Sphere objects. eps=[1,1], equal semi-axis length

Outputs:
pkf: Princle kinematic value
s1_geom: Volume, surface area, and integral of mean curvature of sphere s1
s2_geom: Volume, surface area, and integral of mean curvature of sphere s2
%}

if s1.eps(1)~=1 || s1.eps(2)~=1 || s2.eps(1)~=1 || s2.eps(2)~=1
    error('Input is not sphere, eps');
    return;
elseif isequal(s1.a./s1.a(1), ones(1,3))==false || isequal(s2.a./s2.a(1), ones(1,3))==false
    error('Input is not sphere, semi-axis');
    return;
end

%Volume
s1_geom.volume = s1.GetVolume();
s2_geom.volume = s2.GetVolume();

% Surface area
s1_geom.area = 4*pi*s1.a(1)^2;
s2_geom.area = 4*pi*s2.a(1)^2;

% Integral of mean curvature
s1_geom.int_curvature = 4*pi*s1.a(1);
s2_geom.int_curvature = 4*pi*s2.a(1);

% PKF
pkf = 8 * pi^2 * (s1_geom.volume + s2_geom.volume) +...
    2 * pi * (s1_geom.area * s2_geom.int_curvature +...
    s2_geom.area * s1_geom.int_curvature);

s1_geom;
s2_geom;

end