function [prob, a, x_mink, a_T, x_mink_T] = linearChanceConstraintEnlargedSQ(sub1, sub2, xx, Sigmax, isplot)
if nargin == 4
    isplot = false;
end

    tic;
    prob = 0;

    % Space transformation so that Sigma_x = eye(3)
    xx_T = Sigmax^0.5 \ xx;

    % % Given objective, find x_g(gradient, norm) 
    % Constuct MinkSumClosedForm of s1 and s2
    minkSum = EnlargedMinkSumClosedForm(sub1, sub2);

    option = optimoptions('lsqnonlin',...
                'Algorithm', 'levenberg-marquardt',...
                'display', 'none',...
                'FunctionTolerance', 1e-16,...
                'OptimalityTolerance', 1e-16);

    % Set initial value of optimization
    R1_T = Sigmax^0.5 \ sub1.getMeanOrientation();
    xx_in_s1_T = R1_T' * xx;
    psi_0_T = [atan2( xx_in_s1_T(3), norm(xx_in_s1_T(1:2)) ),...
                atan2( xx_in_s1_T(2), xx_in_s1_T(1) )]

    psi_opt = lsqnonlin(@(psi) func_lsq_tangent(psi, minkSum, xx_T, Sigmax), psi_0_T,...
        [], [], option)
%         Solution gradient in local frame of s1. Notice that gradients
%         used here are all defined in the body-fixed frame of s1
    m_opt = sub1.sq.GetGradientsFromSpherical(psi_opt);
    n_opt = m_opt ./ norm(m_opt);;
    n_opt_minksum = minkSum.sub1.R(:,:,1) * n_opt;
    
    % Get x_mink
    x_mink = minkSum.GetMinkSumFromNormal(n_opt_minksum);

    a = n_opt_minksum;
    
    b = a' * x_mink;

    % Get probability by lcc equation
    prob = 1/2 + 1/2*erf( (b - a'*xx)/sqrt(2*a'*Sigmax*a));
    t = toc;
    
    a_T = Sigmax^0.5 * a;
    x_mink_T = Sigmax^0.5 \ x_mink;

    if isplot
        plotPlane(a, x_mink+sub1.sq.tc, 'm');
    end

end

% Least squares optimization cost for lcc-tangent
function F = func_lsq_tangent(psi, minkSum, xx_T, Sigmax)

m = minkSum.sub1.sq.GetGradientsFromSpherical(psi);
n = m ./ norm(m);

n_minkSum = minkSum.sub1.R(:,:,1) * n;

xMink = minkSum.GetMinkSumFromNormal(n_minkSum);
xMink_T = Sigmax^0.5 \ xMink;

F =xx_T - xMink_T;

end