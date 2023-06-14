function int_curv = poly3d_integral_mean_curvature(pts, k)
% poly3d_integral_mean_curvature computes the integral of mean curvature
% from convex hull
%
%  Inputs:
%    pts     : Nx3 array, each row stores a 3D point
%    k       : Face information, i.e. results from convhull(pts)
%
%  Output:
%    int_curv: Integral of mean curvature for a convex polyhedron
%
%  Author:
%    Sipu Ruan, ruansp@jhu.edu, Aug 2020

% Find edges
tri = triangulation(k, pts);
edges_point = tri.edges();
edges_face = cell2mat( tri.edgeAttachments(edges_point) );

int_curv = 0;
for i = 1:size(edges_point,1)
    % Length for each edge
    edge_length = norm(pts(edges_point(i,1),:) - pts(edges_point(i,2),:));
    
    % Compute surface normals of the faces connected by the edge
    face_normal1 = surface_normal_3d(pts(k(edges_face(i,1),1),:),...
        pts(k(edges_face(i,1),2),:), pts(k(edges_face(i,1),3),:));
    face_normal2 = surface_normal_3d(pts(k(edges_face(i,2),1),:),...
        pts(k(edges_face(i,2),2),:), pts(k(edges_face(i,2),3),:));
    
    if isnan( norm(face_normal1) ) || isnan( norm(face_normal2) )
        continue;
    end
    
    % Compute dihedral angle using two surface normals
    if norm(face_normal1 - face_normal2) < 1e-10
        dihedral_angle = pi;
    elseif norm(face_normal1 + face_normal2) < 1e-10
        dihedral_angle = 0;
    else
        dihedral_angle = pi - abs(acos( dot(face_normal1, face_normal2) ));
    end
    
    % Compute mean curvature
    mean_curv = (pi - dihedral_angle) * edge_length;
    int_curv = int_curv + mean_curv;
end

int_curv = 1/2 * int_curv;
end

function n = surface_normal_3d(pt1, pt2, pt3)
n = cross(pt2-pt1, pt3-pt1);
n = n/norm(n);
end