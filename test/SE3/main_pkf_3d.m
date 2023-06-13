close all; clear; clc;
add_path()
warning('off')

% Main file for computing Principal Kinematic Formula in 3D
%% Define superquadric and ellipsoid objects
N_trial = 10;
N_vtx = (4:20)' * ones(1,2);
I1 = nan(N_trial, size(N_vtx,1));
I2 = nan(N_trial, size(N_vtx,1));

for i = 1:length(N_vtx)
    disp(['Num vertex: ', num2str( prod(N_vtx(i,:)) )])
    
    for j = 1:N_trial
        s1 = SuperQuadrics({2+10*rand(1,3), 0.2+1.8*rand(1,2), [0,0],...
            zeros(3,1), [1,0,0,0], N_vtx(i,:)});
        s2 = SuperQuadrics({2+10*rand(1,3), [1,1], [0,0],...
            10*ones(3,1), [1,0,0,0], N_vtx(i,:)});
        
        %% Compute PKF
        % approximately accurate value
        [pkf_true, s1_geom_true, s2_geom_true] = pkf_3d(s1, s2, 0, 0);
        
        % use defined number of vertices
        [pkf_poly, s1_geom_poly, s2_geom_poly] = pkf_3d(s1, s2, 1, 0);
        
        % use FCL ellipsoid object
        [pkf_poly2, s1_geom_poly2, s2_geom_poly2] = pkf_3d(s1, s2, 1, 1);
        
        %% Accuracy
        I1(j,i) = pkf_poly/pkf_true * 100;
        I2(j,i) = pkf_poly2/pkf_true * 100;
    end
    
    disp(['Averaged accuracy metric (sampling on surface): ',...
        num2str( mean(I1(:,i)) ), '%'])
    disp(['Averaged accuracy metric (FCL ellipsoid object): ',...
        num2str( mean(I2(:,i)) ), '%'])
end

%% Plot metric statistics
figure; hold on;
w = 1.5;
plot(N_vtx.^2, mean(I1,1), 'bo-', 'LineWidth', w)
plot(N_vtx.^2, mean(I2,1), 'r*-', 'LineWidth', w)
legend({'Mesh on surfaces', 'FCL ellipsoid'})

figure; hold on;
w = 1.5;

std1 = std(I1);
std2 = std(I2);
shadedErrorBar(N_vtx.^2, mean(I1,1), std1,...
    'lineProps',{'b-o','markerfacecolor','b'})
shadedErrorBar(N_vtx.^2, mean(I2,1), std2,...
    'lineProps',{'r-*','markerfacecolor','r'})
legend({'Mesh on surfaces', 'FCL ellipsoid'})

%% Plot two objects
figure; hold on; axis equal;
s1.PlotShape('g', 0.7);
s2.PlotShape('b', 0.7);

s2_fcl = get_fcl_ellipsoid(s2);
k_fcl = convhull(s2_fcl);
trisurf(k_fcl, s2_fcl(:,1), s2_fcl(:,2), s2_fcl(:,3), 'FaceAlpha', 0.3)

lightangle(gca,45,30);
lighting gouraud;