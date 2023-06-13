function isCollision = collisionCheck2D(points1, points2)
% Check if any vertex of polygon2 is inside polygon1
isCollision = any(inpolygon(points2(:, 1), points2(:, 2), points1(:, 1), points1(:, 2)));

% Check if any vertex of polygon1 is inside polygon2
isCollision = isCollision || any(inpolygon(points1(:, 1), points1(:, 2), points2(:, 1), points2(:, 2)));

% Check if the edges of the polygons intersect
for i = 1:size(points1, 1)
    p1 = points1(i, :);
    p2 = points1(mod(i, size(points1, 1)) + 1, :);
    for j = 1:size(points2, 1)
        q1 = points2(j, :);
        q2 = points2(mod(j, size(points2, 1)) + 1, :);
        isCollision = isCollision || checkEdgeIntersection(p1, p2, q1, q2);
    end
end

% Function to check if two line segments intersect
function isIntersecting = checkEdgeIntersection(p1, p2, q1, q2)
    v1 = p2 - p1;
    v2 = q2 - q1;
    cross1 = cross([q1-p1, 0], [v1, 0]);
    cross2 = cross([q2-p1, 0], [v1, 0]);
    cross3 = cross([p1-q1, 0], [v2, 0]);
    cross4 = cross([p2-q1, 0], [v2, 0]);
    
    isIntersecting = (dot(cross1, cross2)<= 0) && (dot(cross3, cross4) <= 0);
end
end
