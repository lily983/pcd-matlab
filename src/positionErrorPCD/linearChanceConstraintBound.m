function [prob, t, a, x_mink, r] = linearChanceConstraintBound(s1, s2, mx, Sigma, method, isplot)
% This function reproduce methods using linear chance constraint (LCC)
% to get PCD value for Gaussian distributed variable. The code is inspired
% based on papers "Chance-Constrained Collision Avoidance for MAVs in
% Dynamic Environments" and "Tight Collision Probability for UAV Motion
% Planning in Uncertain Environment" 
% We first proposed to use LCC on superquadraics (SQ) and to use three ways
% to get direction vector n and distance b
    % Inputs:
        % s1, s2: sphere or ellipsoid or superquadraics
        % Sigma: covariance matrix of position error distribution
        % method: ways to get direction vector n and distance b.
            % 'center-point': vector connecting centers of s1, s2. s2.tc - s1.tc
            % 'closed-point': norm at closed point between s1 s2
            % 'tangent-point': closed point from position error center(s2.tc-s1.tc) to transformed
            % minkowski sum surface(transform space so that Sigma is eye(3))
    % Outputs:
        % prob: PCD value
        % t: computation time
if nargin == 3
    method = 'center-point';
    isplot = false;
elseif nargin == 4
    isplot = false;
end

tic;
prob = 0;

% Get direction vector a
%  Get distance b
switch method
    case 'center-point'
        a = (s2.tc-s1.tc)./norm((s2.tc-s1.tc));
        mink = MinkSumClosedForm(s1,s2,quat2rotm(s1.q),quat2rotm(s2.q));
        x_mink = mink.GetMinkSumFromNormal(quat2rotm(s1.q) \ a)+s1.tc;
    case 'closed-point'
        [flag, ~, pt_cls, condition] = collision_cfc(s1, s2, 'constrained');
        if flag && condition<1e-02
          prob = 1;
          t = toc;
          return;
        elseif condition>1e-02
          error('LCC-closed-point collision_cfc_constrained not converge, failed to find closed-points');
        end
        a = (s2.tc - pt_cls.mink) ./ norm((s2.tc - pt_cls.mink));
        x_mink = pt_cls.mink;
    case 'tangent-point'
        [flag, ~, ~, condition] = collision_cfc(s1, s2, 'constrained');
        if flag && condition<1e-02
          prob = 1;
          t = toc;
          return;
        elseif condition>1e-02
          error('LCC-closed-point collision_cfc_constrained not converge, failed to find closed-points');
        end
        % Notice that n is in canonical frame
        n0 = (s2.tc-s1.tc)./norm((s2.tc-s1.tc));
        minkObj = MinkSumClosedForm(s1,s2,Sigma^0.5 \ quat2rotm(s1.q), Sigma^0.5 \ quat2rotm(s2.q));
        option = optimoptions('fmincon', 'Algorithm', 'interior-point',...
            'display', 'iter');
        n_opt = fmincon(@(n) func_con(n, minkObj, Sigma), n0,...
            [], [], [], [], [], [], @(n) nlcon(n, s1), option);
        x_mink = minkObj.GetMinkSumFromNormal(n_opt);
        % x_mink is in rotated and linear transformed space, we need to do
        % inverse linear transformation to convert it to the rotated space
        x_mink = Sigma^0.5 * x_mink +s1.tc;
        % now we need to transform n_opt to rotated frame
        n_final = quat2rotm(s1.q) * n_opt;
        a = n_final./norm(n_final);
    case 'tangent-point-ellip'
        % this only for two ellipsoids
        if strcmp(getObjectType(s1), 'ellip')==0 && strcmp(getObjectType(s1), 'ellip')==0
            error('lcc-tangent-point-ellip only for two ellipsoids cases');
        end
        % Minkowski sum for two ellipsoids and position error
        % ellipsoid(function of distance r)
        n = sym('x', [1,3]);
        syms r
        % get defining matrix of ellip 1 and 2
        A1 = quat2rotm(s1.q) * diag(s1.a) * quat2rotm(s1.q)';
        A2 = quat2rotm(s2.q) * diag(s2.a) * quat2rotm(s2.q)';
        A3 = Sigma;
        eqn1=(A1^2 * n') / norm(A1 * n') +  (A2^2 * n') / norm(A2 * n') +  ((A3.*r)^2 * n') / norm((A3.*r) * n') == s2.tc;
        eqn2 = norm(n) == 1;
        solution = vpasolve([eqn1; eqn2], [n, r]);
        a = double([solution.x1; solution.x2; solution.x3]);
        r = double(solution.r);
        x_mink = (A1^2 *a) / norm(A1 * a) +  (A2^2 * a) / norm(A2 * a) + s1.tc;
end

b = norm(x_mink-s1.tc);

prob = 1/2 + 1/2*erf( (b-a'*mx)/sqrt(2*a'*Sigma*a));

t = toc;

if isplot
    if strcmp(method, 'tangent-point')
        color = 'r';
    elseif strcmp(method, 'closed-point')
        color = 'b';
    elseif strcmp(method, 'tangent-point-ellip')
        color = 'y';
    elseif strcmp(method, 'center-point')
        color = 'g';
    end
%     visualize_bounding_ellip(s1,s2);
    plotPlane(a, x_mink, color);
end

end

function [m_new, flag] = fixed_point(m_init, minkObj)
tol = 1e-10;
p0 = minkObj.s2.tc;
R1 = quat2rotm(minkObj.s1.q);
flag = 1;
%disp('Start fixed_point iteration');
m01 = m_init;

mink = minkObj.GetMinkSumFromGradient(m01) + minkObj.s1.tc;
m02 = minkObj.s1.GetGradientsFromDirection( ...
    m01 + R1'*(p0-mink)/norm(p0-mink) * norm(m01) );

m_fp = [m01, m02];

for i = 1:100
    % calculate the next two guesses for the fixed point.
    m_new = fixed_point_iteration(m01, m02, minkObj);
    m_fp(:,i+1) = m_new;

    % Test convergence
    if norm(m_new-m02)<tol
        flag = 0;
%        disp('Finish fixed-point iteration, find the final gradient')
        return;
    end
    
    % Update gradients
    m01 = m02;
    m02 = m_new;
end
end

% Fixed point iteration step
function m1 = fixed_point_iteration(m01, m02, minkObj)
R1 = quat2rotm(minkObj.s1.q);

mink1 = R1'*minkObj.GetMinkSumFromGradient(m01);
mink2 = R1'*minkObj.GetMinkSumFromGradient(m02);

n01 = m01/norm(m01);
n02 = m02/norm(m02);
r1 = -(n02 - (n02'*n01)*n01)' * (mink2-mink1) /...
    (cross(n02, n01)' * cross(n02, n01)) * n02;

m1 = minkObj.s1.GetGradientsFromDirection(...
    R1'*(minkObj.s2.tc - minkObj.s1.tc) - (mink2 + r1) );
end

function plotPlane(a, x_mink, color)
[x, y]= meshgrid(-1.5:0.1:1.5);

zPlane = -a(1)/a(3) * x - a(2)/a(3)*y + a(1)/a(3)*x_mink(1) + a(2)/a(3)*x_mink(2) + x_mink(3);

scatter3(x_mink(1), x_mink(2), x_mink(3), 'MarkerFaceColor', color);

surf(x, y, zPlane, 'FaceColor', color, 'FaceAlpha', 0.5);
end

% Distance cost: distance from point to Minkowski sums boundary
function F = func_con(n, minkObj, Sigma)
% Minkowski sums with parameter x
p0 = Sigma^0.5 \ minkObj.s2.tc;
mink = minkObj.GetMinkSumFromNormal(n) +  Sigma^0.5 \ minkObj.s1.tc;

% Cost function
F = sum((p0 - mink).^2);
end

% Nonlinear constraint for gradient vector
function [c,ceq] = nlcon(n, s1)
x =  s1.GetPointsFromNormal(n);
ceq = s1.GetImplicitFunction(x);
c = [];
end