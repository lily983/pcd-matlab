function e = metric_fitting_multi_Strue(S_query, S_true_cells, res)
    % IoU between a query convex body and the union of multiple true bodies.
    %
    % INPUTS:
    %   S_query      - 3xN matrix (query body points)
    %   S_true_cells - cell array of {3xMi} (true bodies' points)
    %   res          - grid resolution (optional, default: 20)
    %
    % OUTPUT:
    %   e - IoU = Vol( S_query ∩ (⋃ S_true_i) ) / Vol( S_query ∪ (⋃ S_true_i) )

    if nargin < 3, res = 20; end

    % --- Input checks ---
    if isempty(S_query) || size(S_query,1) ~= 3
        error('S_query must be a non-empty 3xN matrix.');
    end
    if ~iscell(S_true_cells)
        error('S_true_cells must be a cell array of 3xMi matrices.');
    end
    for i = 1:numel(S_true_cells)
        Pi = S_true_cells{i};
        if ~isempty(Pi) && size(Pi,1) ~= 3
            error('Each entry of S_true_cells must be a 3xMi matrix or empty.');
        end
    end

    % --- Aggregate points for bbox ---
    all_true_points = [];
    for i = 1:numel(S_true_cells)
        Pi = S_true_cells{i};
        if ~isempty(Pi)
            all_true_points = [all_true_points, Pi]; %#ok<AGROW>
        end
    end
    if isempty(all_true_points)
        warning('S_true_cells union is empty. Returning IoU = 0.');
        e = 0; return;
    end

    all_points = [S_query, all_true_points];

    % --- Padded bounding box (scale with scene size) ---
    min_coords = min(all_points, [], 2);
    max_coords = max(all_points, [], 2);
    extent = max_coords - min_coords;
    pad = 0.05 * max(extent);                  % 5% of max dimension
    if pad == 0, pad = 1e-2; end               % fallback tiny pad
    min_coords = min_coords - pad;
    max_coords = max_coords + pad;

    % --- Uniform grid over bbox ---
    [x, y, z] = ndgrid( ...
        linspace(min_coords(1), max_coords(1), res), ...
        linspace(min_coords(2), max_coords(2), res), ...
        linspace(min_coords(3), max_coords(3), res));
    grid_points = [x(:)'; y(:)'; z(:)'];

    % --- Membership: query ---
    inside_query = inside_from_points(S_query, grid_points);

    % --- Membership: union of true bodies ---
    inside_true_union = false(size(grid_points,2),1);
    for i = 1:numel(S_true_cells)
        Pi = S_true_cells{i};
        if isempty(Pi), continue; end
        inside_true_union = inside_true_union | inside_from_points(Pi, grid_points);
    end

    % --- IoU via voxel counts ---
    inside_intersection = inside_query & inside_true_union;
    inside_union        = inside_query | inside_true_union;

    denom = sum(inside_union);
    if denom == 0
        % If we reach this, all membership tests reported empty interiors.
        % Most likely inputs are coplanar/degenerate for 3-D.
        warning('Empty union region detected after robustness passes. Returning IoU = 0.');
        e = 0; return;
    end

    e = sum(inside_intersection) / denom;
end

% ---------- Helpers ----------

function inside = inside_from_points(P, grid_points)
    % Robust 3-D "inside" classification from a point cloud approximating a solid.
    % Tries in order:
    %   1) delaunayTriangulation + pointLocation
    %   2) alphaShape with alpha=Inf (convex hull membership)
    %   3) thin-thickness buffer heuristic for near-coplanar sets

    % Deduplicate to avoid unstable triangulations
    P = unique(P','rows')';
    inside = false(size(grid_points,2),1);

    % Quick reject
    if size(P,2) < 4
        % Fewer than 4 points cannot span a volume in 3-D
        inside = false(size(grid_points,2),1);
        return;
    end

    % Check rank (degeneracy)
    Pc = P - mean(P,2);
    r = rank(Pc, 1e-9);

    % --- Try 1: Delaunay 3-D tessellation if full-rank ---
    if r >= 3
        try
            DT = delaunayTriangulation(P');
            [ID, ~] = DT.pointLocation(grid_points');
            inside = ~isnan(ID);
            if any(inside), return; end
        catch
            % fall through
        end
    end

    % --- Try 2: alphaShape convex hull membership ---
    try
        shp = alphaShape(P(1,:)', P(2,:)', P(3,:)', Inf); % convex hull
        inside = inShape(shp, grid_points(1,:)', grid_points(2,:)', grid_points(3,:)');
        if any(inside), return; end
    catch
        % fall through
    end

    % --- Try 3: thin-thickness buffer heuristic (for near-coplanar sets)
    % Project to best-fit plane, define small thickness around it.
    try
        [U,~,~] = svd(Pc, 'econ');
        n = U(:,3);                           % approx normal
        thickness = 1e-3 * max(range(P,2));   % tiny buffer relative to size
        d = abs(n' * (grid_points - mean(P,2)));
        inside = d < thickness;               % classify points close to the sheet as inside
    catch
        inside = false(size(grid_points,2),1);
    end
end
