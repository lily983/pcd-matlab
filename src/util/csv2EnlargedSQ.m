function csvstring = csv2EnlargedSQ(csv, distribution, CL)
%CL: confidence level('95', '99')

if nargin <3
    CL='95';
end

%axis-angle
SQaxis = csv(1, 9:11);
SQaxis = SQaxis./norm(axis,2);
SQangle = csv(1, 12);
quat = axang2quat([SQaxis, SQangle]);

% Build original SQ
sq = SuperQuadrics({csv(1, 1:3), [csv(4), csv(5)], [0,0],...
    [csv(6), csv(7), csv(8)]', [quat(1), quat(2), quat(3), quat(4)], [10, 10]});

% enlarged SQ
Mu = [eye(3), sq.tc; 0, 0, 0, 1];
Sigma = zeros(6,6);
Sigma(4:6, 4:6) = eye(3,3);

% position error covariance
if distribution=="large"
    Sigma(4,4)=1.6e-03;
    Sigma(5,5)=1.4e-03;
    Sigma(6,6)=1.8e-03;
elseif distribution=="small"
    Sigma(4,4)=2.1e-04;
    Sigma(5,5)=2.4e-04;
    Sigma(6,6)=2.7e-04;
end

enlarged_SQ = enlarged_bounding_superquadrics(sq, Mu, Sigma, CL);

enlarged_SQ_csv(1, 1:3) = enlarged_SQ.a;
enlarged_SQ_csv(1, 4:5) = enlarged_SQ.eps;
enlarged_SQ_csv(1, 6:8) = enlarged_SQ.tc';

% quat
enlarged_SQ_axang = quat2axang(enlarged_SQ.q);
enlarged_SQ_csv(1, 9:11) = enlarged_SQ_axang(1,1:3);
enlarged_SQ_csv(1, 12) = enlarged_SQ_axang(1,4);

% visualization
figure; hold on; axis equal;
sq.PlotShape('r', 0.6);
enlarged_SQ.PlotShape('b', 0.3);
visualize_position_error(sq, sq, Sigma(4:6, 4:6));

%print out csv file
csvstring="";
for i=1:size(enlarged_SQ_csv,2)
    csvstring=csvstring+string(enlarged_SQ_csv(i));
    csvstring = csvstring+",";
end

end