function sq=csv2SQ(csv)

%axis-angle
SQaxis = csv(1, 9:11);
SQaxis = SQaxis./norm(axis,2);
SQangle = csv(1, 12);
quat = axang2quat([SQaxis, SQangle]);

% Build original SQ
sq = SuperQuadrics({csv(1, 1:3), [csv(4), csv(5)], [0,0],...
    [csv(6), csv(7), csv(8)]', [quat(1), quat(2), quat(3), quat(4)], [20, 20]});

Sigma = zeros(6,6);
Sigma(4,4)=4.8e-04;
Sigma(5,5)=4.8e-04;
Sigma(6,6)=6e-04;

% Sigma= [0.2546, -0.1820, 0.0208;
%             -0.1820, 0.3172, 0.1886;
%             0.0208, 0.1886, 0.6166];

% visualization
% figure; hold on; axis equal;
% sq.PlotShape('r', 0.6);
visualize_position_error(sq, sq, Sigma(4:6, 4:6));
% visualize_position_error(sq, sq, Sigma);

end