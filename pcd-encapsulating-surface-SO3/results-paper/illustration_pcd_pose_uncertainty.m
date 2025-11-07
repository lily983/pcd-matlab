
clc;clear;close all

%% Loading rotations from samples
T1 = csv2group(tinychair, 'PCG3');
R1 = T1(1:3,1:3,:);

T2 = csv2group(woodenbowl1, 'PCG3');
R2 = T2(1:3,1:3,:);
%% Loading obj parameters
% axes = [0.034,0.021, 0.016];
axes1=[0.27727629926627784,0.2358111973954244,0.17633288963473515]'; 
%% wooden bowl
axes2 = [0.2927352017397958,0.29158779787962896,0.12707043772973667]';
%%
S1  = SuperQuadrics({axes1./2, [0.1,1.0], [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

S2  = SuperQuadrics({axes2./2, [0.2,0.2], [0,0], zeros(3,1), rotm2quat(eye(3)), [20,20]});

figure; hold on
S1.PlotShape('b', 0.1,0.1);
S2.PlotShape('g', 0.1,0.1);

%% 
figure; hold on; axis equal

view([62.1, 32.2822]);    % Azimuth, Elevation
camproj('orthographic');  % Projection type
axis vis3d;  

for i=1:size(T1,3)
    obj_1i = SuperQuadrics({S1.a, S1.eps, [0,0], T1(1:3,4), rotm2quat(R1(:,:,i)), [20,20]});
    obj_1i.PlotShape('b', 0.1,0.2);
    
    obj_2i = SuperQuadrics({S2.a, S2.eps, [0,0], T2(1:3,4), rotm2quat(R2(:,:,i)), [20,20]});
    obj_2i.PlotShape('g', 0.1,0.2);
    
    drawnow;
    pause(0.5);

    if mod(i,1) == 0
        cla; % clear every 2 plots
    end
end
