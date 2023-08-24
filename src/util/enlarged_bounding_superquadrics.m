function s3 = enlarged_bounding_superquadrics(s2, Mu, Sigma, confidenceLevel)
% enlarged_superquadrics_contact_probability This function computes an
% enlarged bounding convex hull of s2 under pose error
% Every variables in the vee vector in se(3) range in mu+/-3*sqrt(sigma)
% when confidence interval is 99% (3*sigma), 95% (2*sigma)

s3 = SuperQuadrics({s2.a, s2.eps, s2.taper, s2.tc, s2.q, s2.N});

mu = get_vee_vector(Mu);

switch confidenceLevel
    case '99'
        CL=3.0;
    case '95'
        CL=2.0;
end

% Get devication 
k=1;
deviation = zeros(6,1);
for i=1:6
    for j=1:6
        if i== j
            deviation(k) = sqrt(Sigma(i,i))*CL;
            k = k+1;
        end
    end
end

vee = zeros(6,64);

vee(:,1) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,2) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,3) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,4) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,5) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,6) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,7) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,8) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,9) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,10) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,11) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,12) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,13) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,14) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,15) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,16) = [mu(1)+deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,17) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,18) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,19) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,20) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,21) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,22) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,23) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,24) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,25) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,26) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,27) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,28) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,29) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,30) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,31) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,32) = [mu(1)+deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,33) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,34) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,35) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,36) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,37) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,38) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,39) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,40) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,41) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,42) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,43) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,44) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,45) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,46) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,47) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,48) = [mu(1)-deviation(1); mu(2)+deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,49) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,50) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,51) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,52) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,53) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,54) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,55) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,56) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)+deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,57) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,58) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,59) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,60) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)+deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

vee(:,61) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)+deviation(6)];

vee(:,62) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)+deviation(5); mu(6)-deviation(6)];

vee(:,63) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)+deviation(6)];

vee(:,64) = [mu(1)-deviation(1); mu(2)-deviation(2); mu(3)-deviation(3); mu(4)-deviation(4); mu(5)-deviation(5); mu(6)-deviation(6)];

% Get surface points of s2 at each 64 extreme configurations
points = zeros(3, s2.N(1)*s2.N(2)*64);

% figure; hold on;
% s3.PlotShape('r', 0.2);

for i=1:64
    
    g_matrix = get_SE3_matrix(vee(:,i));
    s3.tc = g_matrix(1:3,4);
    %if norm is less than threshold, then it is position error case
    if norm(deviation(1:3,1)-zeros(3,1))>1e-04
        s3.q = rotm2quat(g_matrix(1:3,1:3));
    end
%     s3.PlotShape('r',0.05)
%     pause(0.2);
    
    start = 1+(i-1)*s2.N(1)*s2.N(2);
    last = s2.N(1)*s2.N(2)+(i-1)*s2.N(1)*s2.N(2);
    
    points(:,start:last) = s3.GetPoints();
end

% scatter3(points(1,:)',points(2,:)', points(3,:)');

% var_init(1:3) = s2.a;
% var_init(4:5) = s2.eps;
% var_init(6:8) = s2.tc';
% var_init(9:12) = s2.q;
% [sq, ~] = superquadric_fitting(points(:,:),var_init);
% bounding_s2 = SuperQuadrics({sq.a, sq.eps, [0,0], sq.pose.tc, sq.pose.q, s2.N});
% bounding_s2.PlotShape('y',0.2)

%convex hull
k = convhull(points');
% sq_patch = trisurf(k,points(1,:)',points(2,:)',points(3,:)','Visible','on', 'FaceAlpha', 0.5); 
sq_patch = trisurf(k,points(1,:)',points(2,:)',points(3,:)','Visible','off');

surf_points1(:,:) = [sq_patch.XData(1,:)', sq_patch.YData(1,:)', sq_patch.ZData(1,:)'];
surf_points2(:,:) = [sq_patch.XData(2,:)', sq_patch.YData(2,:)', sq_patch.ZData(2,:)'];
surf_points3(:,:) = [sq_patch.XData(3,:)', sq_patch.YData(3,:)', sq_patch.ZData(3,:)'];
surf_points(:,:) = [surf_points1(:,:); surf_points2(:,:); surf_points3(:,:)];


% Test if points are on the surface of convex hull
% surf_points = inpolyhedron(sq_patch.Faces(), sq_patch.Vertices(), points(:,:)', 'TOL', 1e-07)

[sq, ~] = superquadric_fitting(surf_points(:,:)');

s3 = SuperQuadrics({sq.a, sq.eps, [0, 0]...
    sq.pose.tc, sq.pose.q, s2.N});

% hold on;
% s2.PlotShape('g',0.4);
% pause(1)
% s3.PlotShape('b', 0.4);

end