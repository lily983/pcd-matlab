function r_def = metric_deficit(S_query, S_true, res)
    % METRIC_DEFICIT Compute modified deficit ratio between two point sets
    %
    % INPUTS:
    %   S_query - 3xN matrix of points from the query convex body
    %   S_true  - 3xM matrix of points from the true convex body
    %   res     - resolution parameter for the volumetric grid (optional, default: 50)
    %
    % OUTPUT:
    %   r_def - modified deficit ratio
    %           r_def = 1 if S_query is a valid upper bound (S_true ⊆ S_query)
    %           r_def < 1 otherwise, decreasing as the under-approximation worsens
    
    if nargin < 3
        res = 20;
    end

    % Check input dimensions
    if size(S_query, 1) ~= 3 || size(S_true, 1) ~= 3
        error('Point sets must be 3xN and 3xM matrices');
    end

    % Get bounding box
    all_points = [S_query, S_true];
    min_coords = min(all_points, [], 2);
    max_coords = max(all_points, [], 2);

    % Add small margin
    margin = 2;
    min_coords = min_coords - margin;
    max_coords = max_coords + margin;

    % Generate uniform grid
    [x, y, z] = ndgrid(linspace(min_coords(1), max_coords(1), res), ...
                       linspace(min_coords(2), max_coords(2), res), ...
                       linspace(min_coords(3), max_coords(3), res));
    grid_points = [x(:)'; y(:)'; z(:)'];

    % Build convex hulls
    DT_query = delaunayTriangulation(S_query');
    DT_true  = delaunayTriangulation(S_true');

    % Check membership
    [ID_query, ~] = DT_query.pointLocation(grid_points');
    [ID_true, ~]  = DT_true.pointLocation(grid_points');

    inside_query = ~isnan(ID_query);
    inside_true  = ~isnan(ID_true);

    % True volume (approximated by number of grid points)
    V_true = sum(inside_true);

    % Deficit volume (true but not in query)
    V_deficit = sum(inside_true & ~inside_query);

    % Modified deficit ratio
    r_def = 1 - (V_deficit / V_true);

end
