function e = metric_fitting(S_query, S_true, res)
    % METRIC_FITTING Compute fitness metric between two point sets
    %
    % INPUTS:
    %   S_query - 3xN matrix of points from the query convex body
    %   S_true  - 3xM matrix of points from the true convex body
    %   res     - resolution parameter for the volumetric grid (optional, default: 50)
    %
    % OUTPUT:
    %   e - fitness metric (ratio of intersection volume to union volume)
    
    % Set default resolution if not provided
    if nargin < 3
        res = 20;
    end
    
    % Check input dimensions
    if size(S_query, 1) ~= 3 || size(S_true, 1) ~= 3
        error('Point sets must be 3xN and 3xM matrices');
    end
    
    % Get bounding box for both shapes
    all_points = [S_query, S_true];
    min_coords = min(all_points, [], 2);
    max_coords = max(all_points, [], 2);
    
    % Add small margin
%     margin = 0.05 * (max_coords - min_coords);
    margin = 2;
    min_coords = min_coords - margin;
    max_coords = max_coords + margin;
    
    % Generate a uniform grid inside the bounding box
    [x, y, z] = ndgrid(linspace(min_coords(1), max_coords(1), res), ...
                      linspace(min_coords(2), max_coords(2), res), ...
                      linspace(min_coords(3), max_coords(3), res));
    
    % Reshape for testing
    grid_points = [x(:)'; y(:)'; z(:)'];
    
%     % Debug: visualize objects
%     figure; hold on; 
%     scatter3(S_query(1,:),S_query(2,:),S_query(3,:));
%     scatter3(S_true(1,:),S_true(2,:),S_true(3,:));
%     scatter3(x(:), y(:), z(:));
    
    % Use MATLAB's built-in functions for checking if points are in convex hull
    % Get DT object for each point set
    DT_query = delaunayTriangulation(S_query');
    DT_true = delaunayTriangulation(S_true');

%         %Debug: visualize if points are in delaunayTriangulation
%         figure; hold on ; scatter3(S_true(1,:),S_true(2,:),S_true(3,:)); tetramesh(DT_true,'FaceAlpha',0.3);
%         figure; hold on ; scatter3(S_query(1,:),S_query(2,:),S_query(3,:)); tetramesh(DT_query,'FaceAlpha',0.3);

    % Check which grid points are inside each shape
    [ID_query, ~] = DT_query.pointLocation(grid_points');
    [ID_true, ~] = DT_true.pointLocation(grid_points');

    ID_intersection = find(~isnan(ID_query) & ~isnan(ID_true));
    ID_union = find(~isnan(ID_query) | ~isnan(ID_true));

%         % Debug: check if the insides points are correct
%         intersection_pts = grid_points(:, ID_intersection);
%         union_pts = grid_points(:, ID_union);
%         figure; hold on; tetramesh(DT_true,'FaceAlpha',0.1);tetramesh(DT_query,'FaceAlpha',0.1); 
%         scatter3(intersection_pts(1,:),intersection_pts(2,:),intersection_pts(3,:))
%         scatter3(union_pts(1,:),union_pts(2,:),union_pts(3,:))

    % Compute fitness metric
    e = size(ID_intersection, 1) / size(ID_union, 1);
   
end