function ellip = csvSQ2BoundingEllip(csv, plotEllip)
if nargin == 1
    plotEllip = 0;
end

ellip = zeros(size(csv));

for i = 1:size(csv, 1)
    %axis-angle
    SQaxis = csv(i, 9:11);
    SQaxis = SQaxis./norm(axis,2);
    SQangle = csv(i, 12);
    quat = axang2quat([SQaxis, SQangle]);

    % Build original SQ
    s1 = SuperQuadrics({csv(i, 1:3), [csv(i,4), csv(i,5)], [0,0],...
        [csv(i,6), csv(i,7), csv(i,8)]', [quat(1), quat(2), quat(3), quat(4)], [10, 10]});
%     s1.PlotShape('b', 0.8);

    % Get outer bounding ellipsoid for SQ
    s1Points = s1.GetPoints();
    [invSigma1, e1] = MinVolEllipse(s1Points, 0.001);
    Sigma1 = inv(invSigma1);
    [R, lamda_power2, ~] = svd(Sigma1);
    axang_ellip = rotm2axang(R);
    
    ellip(i,:) = [sqrt(diag(lamda_power2))', 1, 1, e1', axang_ellip(1, 1:3), axang_ellip(1, 4)];
    
    if plotEllip
        % plot bounding ellips 
        figure; hold on
%         plot_ellipse(Sigma1, e1, 'edgecolor', 'g');
        sq_ellip = SuperQuadrics({sqrt(diag(lamda_power2))', [1, 1], [0,0], e1, rotm2quat(R), [20,20]});
        sq_ellip.PlotShape('y',0.9);
        s1.PlotShape('b', 0.8);
        
%         % plot position error
%         Sigma = positionErrorDistribution('small');
%         
%         x1_rand = mvnrnd(zeros(3,1), Sigma, 20);
%         s3 = SuperQuadrics({s1.a, s1.eps, s1.taper, s1.tc, s1.q, s1.N});
%          figure; hold on;axis equal;axis off
%          s1.PlotShape('g', 0.6);
%          for j = 1:size(x1_rand, 1)
%             s3.tc = s1.tc + x1_rand(j,:)';
%             s3.PlotShape('g', 0.1);
%          end
%          
%         % plot position error for ellip
%         x4_rand = mvnrnd(zeros(3,1), Sigma, 20);
%         s4 = SuperQuadrics({sqrt(diag(lamda_power2))', [1, 1], [0,0], e1, rotm2quat(R), [20,20]});
%          figure; hold on;axis equal;axis off
%          sq_ellip.PlotShape('b', 0.6);
%          for j = 1:size(x4_rand, 1)
%             s4.tc = sq_ellip.tc + x4_rand(j,:)';
%             s4.PlotShape('b', 0.1);
%          end
    end

end

% get csv file for bounding ellipsoids
ellip_table = table(ellip);
writetable(ellip_table, 'bounding_ellip_for_sq.csv');

return