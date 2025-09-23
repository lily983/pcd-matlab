function r_def = metric_deficit_multi_Strue(S_query, S_true_cells, res)
    % Modified deficit ratio for a query body vs. the UNION of multiple true bodies.
    % r_def = 1 if S_true ⊆ S_query; r_def < 1 otherwise.

    if nargin < 3, res = 20; end

    % --- Input checks ---
    if isempty(S_query) || size(S_query,1) ~= 3
        error('S_query must be a non-empty 3xN matrix.');
    end
    if ~iscell(S_true_cells)
        error('S_true_cells must be a cell array of 3xMi matrices.');
    end
    % Dedup & basic shape checks
    S_query = unique(S_query','rows')';
    for i = 1:numel(S_true_cells)
        Pi = S_true_cells{i};
        if isempty(Pi), S_true_cells{i} = []; continue; end
        if size(Pi,1) ~= 3, error('Each true body must be 3xMi.'); end
        S_true_cells{i} = unique(Pi','rows')';
    end

    % --- Aggregate points for bbox ---
    all_true_points = [];
    for i = 1:numel(S_true_cells)
        if ~isempty(S_true_cells{i})
            all_true_points = [all_true_points, S_true_cells{i}]; %#ok<AGROW>
        end
    end

    % If there are literally no true points, define r_def = 1 (empty set ⊆ S_query)
    if isempty(all_true_points)
        warning('metric_deficit_multi_Strue:EmptyTrueSet', ...
                'True set is empty. Returning r_{def} = 1 by definition.');
        r_def = 1; 
        return;
    end

    all_points = [S_query, all_true_points];
    min_coords = min(all_points, [], 2);
    max_coords = max(all_points, [], 2);
    extent = max_coords - min_coords;
    pad = 0.05 * max(extent); 
    if pad == 0, pad = 1e-2; end
    min_coords = min_coords - pad;
    max_coords = max_coords + pad;

    % --- Grid ---
    [x, y, z] = ndgrid( ...
        linspace(min_coords(1), max_coords(1), res), ...
        linspace(min_coords(2), max_coords(2), res), ...
        linspace(min_coords(3), max_coords(3), res));
    grid_points = [x(:)'; y(:)'; z(:)'];

    % --- Membership: query ---
    inside_query = inside_from_points(S_query, grid_points);

    % --- Membership: union of true bodies ---
    inside_true = false(size(grid_points,2),1);
    for i = 1:numel(S_true_cells)
        Pi = S_true_cells{i};
        if isempty(Pi), continue; end
        inside_true = inside_true | inside_from_points(Pi, grid_points);
    end

    % --- Voxels/volumes ---
    V_true    = sum(inside_true);
    V_deficit = sum( inside_true & ~inside_query );

    % Handle empty true volume (prevents 0/0 -> NaN)
    if V_true == 0
        % No true voxels recognized -> treat as empty union; return 1 by definition.
        warning('metric_deficit_multi_Strue:ZeroTrueVolume', ...
                'No true voxels detected (likely degenerate/copla\-nar). Returning r_{def} = 1.');
        r_def = 1;
        return;
    end

    % Modified deficit ratio
    r_def = 1 - (V_deficit / V_true);
end

% ---------- Helper (robust membership) ----------
function inside = inside_from_points(P, grid_points)
    % Robust 3D "inside" classifier from a point set approximating a solid.
    % Order: Delaunay (if full-rank) -> alphaShape(alpha=Inf) -> thin-sheet fallback.
    inside = false(size(grid_points,2),1);
    if size(P,2) < 4, return; end

    Pc = P - mean(P,2);
    rnk = rank(Pc, 1e-9);

    % Try 1: 3D tetrahedralization
    if rnk >= 3
        try
            DT = delaunayTriangulation(P');
            [ID, ~] = DT.pointLocation(grid_points');
            inside = ~isnan(ID);
            if any(inside), return; end
        catch
            % fall through
        end
    end

    % Try 2: convex hull via alphaShape
    try
        shp = alphaShape(P(1,:)', P(2,:)', P(3,:)', Inf); % convex hull
        inside = inShape(shp, grid_points(1,:)', grid_points(2,:)', grid_points(3,:)');
        if any(inside), return; end
    catch
        % fall through
    end

    % Try 3: thin-sheet fallback (nearly coplanar point sets)
    try
        [U,~,~] = svd(Pc, 'econ');
        n = U(:,3); % approx normal
        thickness = 1e-3 * max(range(P,2)); % tiny buffer
        d = abs(n' * (grid_points - mean(P,2)));
        inside = d < thickness;
    catch
        inside = false(size(grid_points,2),1);
    end
end
