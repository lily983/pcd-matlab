function visualize_bounding_ellip(s1, s2)
% Plot the minkowski sum of s1 and s2 and the bounding ellip of the
% minkowski sum

[mf, Sigmaf] = get_bounding_ellip(s1, s2);
% 
% % [mf, Sigmaf] = get_bounding_ellip_fixed_point(s1, s2);

% invT = inv(Sigmaf)+inv(Sigmax);
% % M = 1/det(T).*T;
% c = (det(Sigmaf)/det(inv(invT)))^(1/3);
% T = c.*inv(invT);

% R diag R'
[R, Lamda ,~] = svd(Sigmaf);
    
dimension = size(s1.a,2);
    
    if dimension == 2
        s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
            s2.tc(1), s2.tc(2), s2.ang, s2.N]);
        
        %s3 is the bounding ellip
        s3.a = sqrt(diag(Lamda));
        s3.ang = rotm2angle(R);
        s3.tc = mf;
        
        % Get minkowski sum of s1 and s2
        s1Points = s1.GetPoints()';
        s2Points = -1*s2.GetPoints()';
        pgon1 = polyshape(s1Points(:,1), s1Points(:,2));
        pgon2= polyshape(s2Points(:,1), s2Points(:,2));
        mSum = minkowskiSum(pgon1, pgon2);
        
%         figure; hold on;axis equal; axis off;
        patch(mSum.Vertices(:,1), mSum.Vertices(:,2), hex2rgb('CE2836'), 'FaceAlpha', 0.2, 'EdgeColor',...
            hex2rgb('CE2836'), 'LineWidth', 2, 'LineStyle', '--');
%         s1.PlotShape(hex2rgb('4FAAD1'), 0.6);
%         s2.PlotShape(hex2rgb('45AC59'), 0.6);
%         s3.PlotShape(hex2rgb('EBBF00 '), 0.7);

    elseif dimension ==3
        s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
            s2.tc, s2.q, s2.N});
        s3.a = sqrt(diag(Lamda));
        s3.q = rotm2quat(R);
%         s3.tc = mf;
        m1 = s1.GetGradientsCanonical();
        mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q), quat2rotm(s2.q));
%         x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc-s2.tc;
        x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;
%         mu_m = T* inv(Sigmaf) * mf ;
        s3.tc = s1.tc;

        SN = s1.N;
        
        figure; 
        hold on;
%         axis equal;
        s1.PlotShape(hex2rgb('4FAAD1'), 0.7);
        s2.PlotShape(hex2rgb('45AC59'), 0.7);
        s3.PlotShape(hex2rgb('EBBF00 '), 0.4);
        patch_mink = surf2patch(reshape(x_mink(1,:), SN), reshape(x_mink(2, :), SN), reshape(x_mink(3, :), SN), 'triangles');
        patch_mink.FaceColor = hex2rgb('CE2836 ');
        patch_mink.FaceAlpha = 0.01;
        patch_mink.EdgeColor = hex2rgb('CE2836');
        patch_mink.LineWidth = 1;
        patch_mink.LineStyle = '--';
        patch(patch_mink);
        
    end
     
    
end