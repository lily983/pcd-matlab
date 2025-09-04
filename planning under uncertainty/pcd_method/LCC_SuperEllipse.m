function [prob, t, a, x_mink, m_current, x_mink_T] = LCC_SuperEllipse(s1, s2, xx, Sigmax, method, isplot)
if nargin == 5
    isplot = false;
end

switch method
    case 'center-point'
        
        % Get bounding ellipsoid for each object
        s1Points = s1.GetPoints()-s1.tc;
        s2Points = s2.GetPoints()-s2.tc;
        
        [invSigma1, c1] = MinVolEllipse(s1Points, 0.001);
        [invSigma2, c2] = MinVolEllipse(s2Points, 0.001);
        Sigma1 = inv(invSigma1);
        Sigma2 = inv(invSigma2);
        
        % We don't count in computation time of bounding ellipsoid for
        % objects
        tic;
        prob = 0;
        
        % Get bounding ellip for the Minkowski sum of two ellips
        Sigmaf = (1 + sqrt(trace(Sigma2)/trace(Sigma1))) * Sigma1...
            + (1 + sqrt(trace(Sigma1)/trace(Sigma2))) * Sigma2;
        
        % transformed relative position error
        xx_T = Sigmaf^0.5 \ (xx+c2-c1);
        Sigmax_T = Sigmaf^0.5 \ Sigmax / Sigmaf^0.5;
        
        % Get normal vector after space transformation
        a = xx_T./norm(xx_T);
        % The bounding ellip becomes a sphere and thus the radius is one
        b = 1; 
        
        % Get probability by lcc equation
        prob = 1/2 + 1/2*erf( (b-a'*xx_T)/sqrt(2*a'*Sigmax_T*a));
        t = toc;
        
        % Get norm vector, x_mink in the untransformed space
        a = (Sigmaf \ xx) ./ norm(Sigmaf \ xx);
        x_mink = xx./norm(Sigmaf^0.5 \ xx);
    case 'center-point-cfc'
        % Find x_mink: intersection point between the exact minksum
        % boundary and line connecting xx and origin
        tic;
        prob = 0;
        
        % % Given objective, find x_mink(gradient, norm) 
        % Constuct MinkSumClosedForm of s1 and s2
        minkSum = MinkSumClosedForm(s1, s2, angle2rotm(s1.ang), angle2rotm(s2.ang));

        % the intersection point should lie at the same line as origin to xx 
        xx_norm = xx./norm(xx);

        option = optimoptions('lsqnonlin',...
                    'Algorithm', 'levenberg-marquardt',...
                    'display', 'none',...
                    'FunctionTolerance', 1e-8,...
                    'OptimalityTolerance', 1e-8);
        
        % Set initial value of optimization
        R1 = angle2rotm(s1.ang);
        s2_tc_in_s1 = R1' * (s2.tc-s1.tc);
        psi_0 = atan2( s2_tc_in_s1(1), norm(s2_tc_in_s1(2)));
        
        psi_opt = lsqnonlin(@(psi) func_lsq(psi, minkSum, xx_norm), psi_0,...
            [], [], option);
%         Solution gradient in local frame of s1. Notice that gradients
%         used here are all defined in the body-fixed frame of s1
        m_opt = s1.GetGradientFromAngle(psi_opt);
        
        % Get x_mink
        x_mink = minkSum.GetMinkSumFromGradient(m_opt);
        
        % Get m_opt in current space
        %Noted that m_opt is defined in s1's body frame
        m_current =R1' \ m_opt; 
        a = m_current ./ norm(m_current);
        
        % Get b
        b = x_mink' * a;
        
        % Get probability by lcc equation
        prob = 1/2 + 1/2*erf( (b-a'*xx)/sqrt(2*a'*Sigmax*a));
        t = toc;
    case 'tangent-point-cfc'
        % Find x_g': intersection point between the exact minksum
        % boundary and confidence level surface (after space
        % transformation)
        tic;
        prob = 0;
        
        % Space transformation so that Sigma_x = eye(2)
        xx_T = Sigmax^0.5 \ xx;
%  angle2rotm(s1.ang), angle2rotm(s2.ang)
        % % Given objective, find x_g(gradient, norm) 
        % Constuct MinkSumClosedForm of s1 and s2
        minkSum_T = MinkSumClosedForm(s1, s2, Sigmax^0.5 \ angle2rotm(s1.ang), Sigmax^0.5 \ angle2rotm(s2.ang));

        option = optimoptions('lsqnonlin',...
                    'Algorithm', 'levenberg-marquardt',...
                    'display', 'none',...
                    'FunctionTolerance', 1e-16,...
                    'OptimalityTolerance', 1e-16);
        
        % Set initial value of optimization
        R1_T = Sigmax^0.5 \ angle2rotm(s1.ang);
        s2_tc_in_s1_T = R1_T' * (s2.tc-s1.tc);
        psi_0_T = atan2( s2_tc_in_s1_T(1), norm(s2_tc_in_s1_T(2)));
        
        psi_opt = lsqnonlin(@(psi) func_lsq_tangent(psi, minkSum_T, xx_T), psi_0_T,...
            [], [], option);
%         Solution gradient in local frame of s1. Notice that gradients
%         used here are all defined in the body-fixed frame of s1
        m_opt = s1.GetGradientFromAngle(psi_opt);
        
        % Get x_mink
        x_mink_T = minkSum_T.GetMinkSumFromGradient(m_opt);
        
        % Get m_opt in current space
        %Noted that m_opt is defined in s1's body frame
        m_current = (Sigmax^0.5 \ angle2rotm(s1.ang))' \ m_opt; 
        a_T = m_current ./ norm(m_current);
        
        % Get b
        b_T = x_mink_T' * a_T;
        
        % Get probability by lcc equation
        prob = 1/2 + 1/2*erf( (b_T-a_T'*xx_T)/sqrt(2*a_T'*eye(2)*a_T));
        t = toc;
        m = angle2rotm(s1.ang) * m_opt; 
        a = m ./ norm(m);
        x_mink = Sigmax^0.5 * x_mink_T;
%         return
end

if isplot
    if strcmp(method, 'center-point')
        color = 'g';
    elseif strcmp(method, 'center-point-cfc')
        color = 'm';
    end
%     visualize_bounding_ellip(s1,s2);
    color = 'r';
    plotPlane(a, x_mink, color);
end

end

% Least squares optimization cost
function F = func_lsq(psi, minkSum, xx)
m = minkSum.s1.GetGradientFromAngle(psi);

xMink = minkSum.GetMinkSumFromGradient(m);

F = norm(xMink) - xMink' * xx;
end

% Least squares optimization cost for lcc-tangent
function F = func_lsq_tangent(psi, minkSum_T, xx_T)
m = minkSum_T.s1.GetGradientFromAngle(psi);

xMink = minkSum_T.GetMinkSumFromGradient(m);

% F = 0.5 * sum((xx_T - xMink).^2);
F =xx_T - xMink;
% psi
% xx_T
end