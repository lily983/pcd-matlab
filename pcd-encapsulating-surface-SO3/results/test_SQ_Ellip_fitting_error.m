

close all; clear; clc;
add_paths();

N_trial = 10;

disp('****************************************************************')
disp('*** Testing for geometric model fitting approximation errors ***')
disp('****************************************************************')

%% Fitting superquadrics to convex polyhedra
disp('Fitting superquadrics to convex polyhedra...')
N_pts = 100;

vol_poly = nan(1, N_trial);

vol_sq = nan(1, N_trial);
rel_vol_sq_poly = nan(1, N_trial);
error_rel_sq_poly = nan(1, N_trial);
error_pt_surf = nan(1, N_trial);

vol_ellip = nan(1, N_trial);
rel_vol_ellip_poly = nan(1, N_trial);
error_rel_ellip_poly = nan(1, N_trial);
error_pt_surf_ellip = nan(1, N_trial);

ratio_vol_sq_poly = nan(1, N_trial);
ratio_vol_ellip_poly = nan(1, N_trial);

for i = 1:N_trial
    % Generate random convex polyhedra
    point_cloud = [20;15;10] .* (2*rand(3,N_pts) - 1);
    [k, vol_poly(i)] = convhull(point_cloud(1,:), point_cloud(2,:),...
        point_cloud(3,:));
    poly = point_cloud(:,unique(k));
    
%     figure;hold on
%     trisurf(k, point_cloud(1,:), point_cloud(2,:), point_cloud(3,:))
    
    % Fit superquadric model
    [sq_shape, error_pt_surf(i)] = sq_fitting([], poly);
    sq = SuperQuadrics({sq_shape.a, sq_shape.eps, [0,0],...
        sq_shape.pose.tc, sq_shape.pose.q, [50, 50]});
    vol_sq(i) = sq.GetVolume();
    
    poly_canonical = quat2rotm(sq.q)\(poly - sq.tc);
    error_pt_surf(i) = sum( abs(sq.GetImplicitFunction(poly_canonical)) )...
        / size(poly_canonical, 2);
    
    % Fit ellipsoid model
    [ellip_shape, error_pt_surf_ellip(i)] = ellip_fitting([], poly);
    ellip = SuperQuadrics({ellip_shape.a, ellip_shape.eps, [0,0],...
        ellip_shape.pose.tc, ellip_shape.pose.q, [50, 50]});
    vol_ellip(i) = ellip.GetVolume();
    
    poly_canonical_ellip = quat2rotm(ellip.q)\(poly - ellip.tc);
    error_pt_surf_ellip(i) = sum( abs(ellip.GetImplicitFunction(poly_canonical_ellip)) )...
        / size(poly_canonical, 2);
    
    % Relative volume
    error_rel_sq_poly(i) = abs(vol_sq(i)-vol_poly(i))./vol_poly(i);
    rel_vol_sq_poly(i) = vol_sq(i)./vol_poly(i);
    
    error_rel_ellip_poly(i) = abs(vol_ellip(i)-vol_poly(i))./vol_poly(i);
    rel_vol_ellip_poly(i) = vol_ellip(i)./vol_poly(i);
    
%     % intersection volume / union volume
    ratio_vol_sq_poly(i) = metric_fitting(sq.GetPoints(), poly); 
    ratio_vol_ellip_poly(i) = metric_fitting(ellip.GetPoints(), point_cloud); 
    
end

% disp(['Mean error avg point to fitted surface distance: ',...
%     num2str( mean(error_pt_surf) )]);
disp(['Mean error rel volume SQ fit Polyhedra: ',...
    num2str( mean(error_rel_sq_poly) )]);
% disp(['Mean relative volume SQ fit Polyhedra: ',...
%     num2str( mean(rel_vol_sq_poly) )]);
% disp(['Mean error avg point to fitted surface distance: ',...
%     num2str( mean(error_pt_surf_ellip) )]);
disp(['Mean error rel volume ellip fit Polyhedra: ',...
    num2str( mean(error_rel_ellip_poly) )]);
% disp(['Mean relative volume ellip fit Polyhedra: ',...
%     num2str( mean(rel_vol_ellip_poly) )]);

disp(['Volume ratio SQ fit Polyhedra: ',...
    num2str( mean(ratio_vol_sq_poly) )]);
disp(['Volume ratio ellip fit Polyhedra: ',...
    num2str( mean(ratio_vol_ellip_poly) )]);

% Plots
figure; hold on; axis equal; axis off;
% plot3(point_cloud(1,:), point_cloud(2,:), point_cloud(3,:), 'bo');

plot3(poly(1,:), poly(2,:), poly(3,:), 'r*');
trisurf(k, point_cloud(1,:), point_cloud(2,:), point_cloud(3,:),...
    'FaceColor', 'y', 'FaceAlpha', 0.7)

sq.PlotShape('g', 0.3);
ellip.PlotShape('b', 0.3);

lightangle(gca,45,30);
lighting gouraud;

%% Fitting superquadrics and ellipsoids to basic primitives
disp('****************************************************************')
disp('Fitting superquadrics and convex polyhedra to basic primitives...')

% Cube
disp('** Cube **')
vol_cube = nan(1, N_trial);

vol_poly_cube = nan(1, N_trial);

vol_sq_cube = nan(1, N_trial);
error_rel_sq_cube = nan(1, N_trial);

vol_inscribed_ellip_cube = nan(1, N_trial);
error_rel_inscribed_ellip_cube = nan(1, N_trial);

vol_circumscribed_ellip_cube = nan(1, N_trial);
error_rel_circumscribed_ellip_cube = nan(1, N_trial);

for i = 1:N_trial
    % Define cube
    cube_size = [50,36,20] .* rand(1,3) + 1;
    vol_cube(i) = prod(cube_size);
    
    % Polyhedron (exact)
    vol_poly_cube(i) = vol_cube(i);
    
    % SQ approximation
    sq_cube = SuperQuadrics({cube_size./2, [0.2,0.2], [0,0],...
            [0;0;0], [0,1,0,0], [50, 50]});
    vol_sq_cube(i) = sq_cube.GetVolume();
    
    % inscribed ellipsoid approximation
    inscribed_ellip_cube = SuperQuadrics({cube_size./2, [1, 1], [0,0],...
            [0;0;0], [0,1,0,0], [50, 50]});
    vol_inscribed_ellip_cube(i) = inscribed_ellip_cube.GetVolume();
    
    % circumscribed ellipsoid approximation
    circumscribed_ellip_cube = SuperQuadrics({cube_size./2.*sqrt(3), [1, 1], [0,0],...
            [0;0;0], [0,1,0,0], [50, 50]});
    vol_circumscribed_ellip_cube(i) = circumscribed_ellip_cube.GetVolume();
    
    % Relative volume
    error_rel_sq_cube(i) = abs(vol_sq_cube(i)-vol_cube(i))./...
        vol_cube(i);
    error_rel_inscribed_ellip_cube(i) = abs(vol_inscribed_ellip_cube(i)-vol_cube(i))./...
        vol_cube(i);
    error_rel_circumscribed_ellip_cube(i) = abs(vol_circumscribed_ellip_cube(i)-vol_cube(i))./...
        vol_cube(i);
    
    % cube points
    cube_points = sample_cube_surface(cube_size(1), cube_size(2), cube_size(3), 60);
    
     % Fit ellipsoid model
    [ellip_cube, ~] = ellip_fitting([], cube_points);
    ellip_cube = SuperQuadrics({ellip_cube.a, ellip_cube.eps, [0,0],...
        ellip_cube.pose.tc, ellip_cube.pose.q, [50, 50]});
    vol_ellip(i) = ellip_cube.GetVolume();
    
    % intersection volume / union volume
    ratio_vol_sq_cube(i) = metric_fitting(sq_cube.GetPoints(), cube_points); 
    ratio_vol_circumscribed_ellip_cube(i) = metric_fitting(circumscribed_ellip_cube.GetPoints(), cube_points); 
    ratio_vol_inscribed_ellip_cube(i) = metric_fitting(inscribed_ellip_cube.GetPoints(), cube_points); 
    ratio_vol_ellip_cube(i) = metric_fitting(ellip_cube.GetPoints(), cube_points); 
end

disp(['Mean rel volume SQ fit Cube: ',...
    num2str( mean(error_rel_sq_cube) )]);
% disp(['Mean rel volume inscribed ellip fit Cube: ',...
%     num2str( mean(error_rel_inscribed_ellip_cube) )]);
disp(['Mean rel volume circumscribed ellip fit Cube: ',...
    num2str( mean(error_rel_circumscribed_ellip_cube) )]);
disp(['Volume ratio SQ fit Cube: ',...
    num2str( mean(ratio_vol_sq_cube) )]);
disp(['Volume ratio incribed ellip fit Cube: ',...
    num2str( mean(ratio_vol_circumscribed_ellip_cube) )]);
disp(['Volume ratio circumscribed ellip fit Cube: ',...
    num2str( mean(ratio_vol_circumscribed_ellip_cube) )]);
disp(['Volume ratio ellip fit Cube: ',...
    num2str( mean(ratio_vol_ellip_cube) )]);

figure; hold on; axis equal;
sq_cube.PlotShape('g', 0.3);
inscribed_ellip_cube.PlotShape('b', 0.4);
circumscribed_ellip_cube.PlotShape('y', 0.2);
ellip_cube.PlotShape('r', 0.4);

scatter3(cube_points(1,:), cube_points(2,:), cube_points(3,:))

lightangle(gca,45,30);
lighting gouraud;

%%
% Cylinder
disp('** Cylinder **')
N_poly_pts = 60;
vol_cyl = nan(1, N_trial);

vol_sq_cyl = nan(1, N_trial);
% vol_rel_sq_cyl = nan(1, N_trial);
error_rel_sq_cyl = nan(1, N_trial);

vol_inscribed_ellip_cyl = nan(1, N_trial);
% vol_rel_inscribed_ellip_cyl = nan(1, N_trial);
error_rel_inscribed_ellip_cyl = nan(1, N_trial);

vol_circumscribed_ellip_cyl = nan(1, N_trial);
% vol_rel_circumscribed_ellip_cyl = nan(1, N_trial);
error_rel_circumscribed_ellip_cyl = nan(1, N_trial);

for i = 1:N_trial
    % Define cylinder
    cyl_size = [50,36,20].*rand(1,3) + 1;
    vol_cyl(i) = pi * prod(cyl_size);
    
    % SQ approximation
    sq_cyl = SuperQuadrics({[cyl_size(1:2), cyl_size(3)/2], [0.1,1],...
        [0,0], [0;0;0], [0,1,0,0], [50, 50]});
    vol_sq_cyl(i) = sq_cyl.GetVolume();
    
    % inscribed ellipsoid approximation
    inscribed_ellip_cyl = SuperQuadrics({[cyl_size(1:2), cyl_size(3)/2], [1, 1], [0,0],...
            [0;0;0], [0,1,0,0], [50, 50]});
    vol_inscribed_ellip_cyl(i) = inscribed_ellip_cube.GetVolume();
    
    % circumscribed ellipsoid approximation
    circumscribed_ellip_cyl = SuperQuadrics({[cyl_size(1:2).*sqrt(1.5), cyl_size(3)*sqrt(1.5)/2], [1, 1], [0,0],...
            [0;0;0], [0,1,0,0], [50, 50]});
    vol_circumscribed_ellip_cyl(i) = circumscribed_ellip_cube.GetVolume();
    
    % Polyhedral approximation
    th = -pi:2*pi/(N_poly_pts/2-1):pi;
    cross_section = cyl_size(1:2)'.*[cos(th);sin(th)];
    pts_cyl = [[cross_section; -cyl_size(3)/2 *...
        ones(1,size(cross_section,2))], [cross_section;...
        cyl_size(3)/2 * ones(1,size(cross_section,2))]];
     [k, vol_poly_cyl(i)] = convhull(pts_cyl');
     
      % Fit ellipsoid model
    [ellip_cyl, ~] = ellip_fitting([], pts_cyl);
    ellip_cyl = SuperQuadrics({ellip_cyl.a, ellip_cyl.eps, [0,0],...
        ellip_cyl.pose.tc, ellip_cyl.pose.q, [50, 50]});
    
     % Relative volume
    error_rel_sq_cyl(i) = abs(vol_sq_cyl(i)-vol_cyl(i))./...
        vol_cyl(i);
    error_rel_inscribed_ellip_cyl(i) = abs(vol_inscribed_ellip_cyl(i)-vol_cyl(i))./...
        vol_cyl(i);
    error_rel_circumscribed_ellip_cyl(i) = abs(vol_circumscribed_ellip_cyl(i)-vol_cyl(i))./...
        vol_cyl(i);
    
    % intersection volume / union volume
    ratio_vol_sq_cyl(i) = metric_fitting(sq_cyl.GetPoints(), pts_cyl); 
    ratio_vol_inscribed_ellip_cyl(i) = metric_fitting(inscribed_ellip_cyl.GetPoints(), pts_cyl); 
    ratio_vol_circumscribed_ellip_cyl(i) = metric_fitting(circumscribed_ellip_cyl.GetPoints(), pts_cyl); 
    ratio_vol_ellip_cyl(i) = metric_fitting(circumscribed_ellip_cyl.GetPoints(), pts_cyl); 
end

disp(['Mean rel volume SQ fit Cube: ',...
    num2str( mean(error_rel_sq_cyl) )]);
disp(['Mean rel volume inscribed ellip fit Cube: ',...
    num2str( mean(error_rel_inscribed_ellip_cyl) )]);
disp(['Mean rel volume circumscribed ellip fit Cube: ',...
    num2str( mean(error_rel_circumscribed_ellip_cyl) )]);

disp(['Volume ratio SQ fit Cube: ',...
    num2str( mean(ratio_vol_sq_cyl) )]);
disp(['Volume ratio inscribed ellip fit Cube: ',...
    num2str( mean(ratio_vol_inscribed_ellip_cyl) )]);
disp(['Volume ratio circumscribed ellip fit Cube: ',...
    num2str( mean(ratio_vol_circumscribed_ellip_cyl) )]);
disp(['Volume ratio ellip fit Cube: ',...
    num2str( mean(ratio_vol_ellip_cyl) )]);

figure; hold on; axis equal; axis off;
plot3(pts_cyl(1,:), pts_cyl(2,:), pts_cyl(3,:), 'r*');
trisurf(k, pts_cyl(1,:), pts_cyl(2,:), pts_cyl(3,:),...
    'FaceColor', 'cyan', 'FaceAlpha', 0.3)

sq_cyl.PlotShape('g', 0.3);
inscribed_ellip_cyl.PlotShape('b', 0.4);
circumscribed_ellip_cyl.PlotShape('y', 0.3);

lightangle(gca,45,30);
lighting gouraud;
