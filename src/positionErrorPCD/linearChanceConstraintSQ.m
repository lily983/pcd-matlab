function [prob, t, a, x_mink] = linearChanceConstraintSQ(s1, s2, xx, Sigmax, method, isplot)
if nargin == 4
    isplot = false;
end

switch method
    case 'center-point'
        % Get bounding ellipsoid for each object
        s1Points = s1.GetPoints();
        s2Points = s2.GetPoints();
        
        [invSigma1, ~] = MinVolEllipse(s1Points, 0.001);
        [invSigma2, ~] = MinVolEllipse(s2Points, 0.001);
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
        xx_T = Sigmaf^0.5 \ xx;
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
        minkSum = MinkSumClosedForm(s1, s2, quat2rotm(s1.q), quat2rotm(s2.q));

        % the intersection point should lie at the same line as origin to xx 
        xx_norm = xx./norm(xx);

        option = optimoptions('lsqnonlin',...
                    'Algorithm', 'levenberg-marquardt',...
                    'display', 'none',...
                    'FunctionTolerance', 1e-8,...
                    'OptimalityTolerance', 1e-8);
        
        % Set initial value of optimization
        psi_0 = [0,0];
        
        psi_opt = lsqnonlin(@(psi) func_lsq(psi, minkSum, xx_norm), psi_0,...
            [], [], option);
%         Solution gradient in local frame of s1. Notice that gradients
%         used here are all defined in the body-fixed frame of s1
        m_opt = s1.GetGradientsFromSpherical(psi_opt);
        
        % Get x_mink
        x_mink = minkSum.GetMinkSumFromGradient(m_opt);
        
        % Get m_opt in current space
        %Noted that m_opt is defined in s1's body frame
        m_current = quat2rotm(s1.q)' \ m_opt; 
        a = m_current ./ norm(m_current);
        
        % Get b
        b = x_mink' * a;
        
        % Get probability by lcc equation
        prob = 1/2 + 1/2*erf( (b-a'*xx)/sqrt(2*a'*Sigmax*a));
        t = toc;
        return
end

if isplot
    if strcmp(method, 'center-point')
        color = 'g';
    elseif strcmp(method, 'center-point-cfc')
        color = 'm';
    end
%     visualize_bounding_ellip(s1,s2);
    plotPlane(a, x_mink+s1.tc, color);
end

end

% Least squares optimization cost
function F = func_lsq(psi, minkSum, xx)
m = minkSum.s1.GetGradientsFromSpherical(psi);

xMink = minkSum.GetMinkSumFromGradient(m);

F = norm(xMink) - xMink' * xx;
end