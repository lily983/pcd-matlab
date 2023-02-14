function [x, p, cost, sigma2] = SuperquadricsGaussian(point, para)

% Written by Weixiao Liu, PhD student @ JHU, NUS
%            Yuwei Wu @ NUS
% May 24th, 2021, Singapore
% -------------------------------------------------------------------------
% DESCRIPTION: This algorithm solves for the optimal superquadrics fitting of a given
%              point cloud. Probabilistic model is adpot to formulate the problem, and
%              thus is roubust enough to tolerate some amount of outliers. The outlier
%              probability is treated as a hidden random variable and is updated via
%              the Bayes' rule (EM algorithm). The parameters of the superquadric is
%              solved iteratively via maximum likelihood estimation.
%
% INPUT: point - point cloud array (3 x N)
%        para - structure include tunable parameters such as w, iterEM_max,
%               tolrenceEM, relative_toleranceEM
%
% OUTPUT: x - fitted superquadrics parameters
%         p - outlier probability of the corresponding points
%         cost - the final value of the cost function
%         sigma2 - the final covariance of the Gaussian components
% -------------------------------------------------------------------------
%% Configuration

% set parameters
w = para.w;
iterEM_max = para.iterEM;
iterEM_min = 2; % 5
toleranceEM = para.toleranceEM;
relative_toleranceEM = para.relative_toleranceEM;
debug_mode = para.debug_mode;
arclength = para.arclength;
adaptive_upper = para.adaptive_upper;

% optimization settings
options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'Display', 'off', 'MaxIterations', para.iterLSQ_max);

%% Initialization

% translate the coordinate to the center of mass
t0 = mean(point, 2);
point = point - t0;

% rescale 
if para.rescale == 1
   max_length = max(max(point));
   scale = max_length / 10; %10
   point = point / scale;
end

% eigen analysis (principal component analysis) for initializing rotation
[EigenVector, ~] = EigenAnalysis(point);
EigenVector = [EigenVector(:, 1), EigenVector(:, 3), cross(EigenVector(:, 1), EigenVector(:, 3))];
euler0 = rotm2eul(EigenVector);

% initialize as median along transformed axises
point_rot0 = EigenVector' * point;
s0 = [median(abs(point_rot0(1, :))), median(abs(point_rot0(2, :))), median(abs(point_rot0(3, :)))];

% initial configuration
x0 = [1, 1, reshape(s0, [1, 3]), euler0, zeros(1, 3)];

% set lower and upper bounds for the superquadrics
upper = 4 * max(max(abs(point)));
lb = [0.0 0.0 0.001 0.001 0.001 -2*pi -2*pi -2*pi -ones(1, 3) * upper]; %0.01 0.001
ub = [2.0 2.0 ones(1, 3) * upper  2*pi 2*pi 2*pi ones(1, 3) * upper];

%------------------------------debug block---------------------------------
if debug_mode == 1
    figure
    [point_fit] = sphericalProduct_sampling(x0, arclength);
    plot3(point(1, :), point(2, :), point(3, :), '.', 'Color', 'r')
    hold on
    plot3(point_fit(1, :), point_fit(2, :), point_fit(3, :), '.', 'Color', 'b')
    title('Initialization')
    axis equal
    hold off
    disp('Initial configuration is:')
    disp(x0)
    waitforbuttonpress
end
%--------------------------------------------------------------------------

%% EMS Algorithm

% set bounding volume of outlier space
% V = (max(point(1, :)) - min(point(1, :))) * (max(point(2, :)) ...
%     - min(point(2, :))) * (max(point(3, :)) - min(point(3, :)));
V = (max(point_rot0(1, :)) - min(point_rot0(1, :))) * (max(point_rot0(2, :)) ...
    - min(point_rot0(2, :))) * (max(point_rot0(3, :)) - min(point_rot0(3, :)));
% outlier probability density
p0 = 1 / V;

% initialize variance for gaussian model
if para.sigma2 == 0
    sigma2 = V ^ (1 / 3) / 10;
else
    sigma2 = para.sigma2;
end

% initialize parameters
x = x0;
cost = 0;
switched = 0;
p = ones(1, size(point, 2));

for iterEM = 1 : iterEM_max
    % evaluating distance
    [dist] = distance(point, x);
    % evaluating corespondence probability
    if w ~= 0
        p = correspendence(dist, sigma2, w, p0);
    end
    % adaptive upper bound
    if  adaptive_upper == 1
        R_current = eul2rotm(x(6 : 8));
        point_rot_current = R_current' * point - R_current' * x(9 : 11)';
%         point_rot_current = R_current' * point
        ub_a = 1.1 * [max(abs(point_rot_current(1, :))), max(abs(point_rot_current(2, :))), max(abs(point_rot_current(3, :)))];
        ub = [2.0 2.0 ub_a  2*pi 2*pi 2*pi ub_a];
        lb = [0.0 0.0 0.01 0.01 0.01 -2*pi -2*pi -2*pi -ub_a];
    end
    % optimization
    cost_func = @(x) weighted_dist(x, point, p, sigma2);
    [x_n, cost_n] = lsqnonlin(cost_func, x, lb, ub, options);
    % update sigma
    sigma2_n = cost_n / (3 * sum(p));
%     w = 1 - sum(p);
    % evaluate relative cost decrease
    relative_cost = (cost - cost_n) / cost_n;
    
    %------------------------------debug block---------------------------------
    if debug_mode == 1
        figure
        [point_fit] = sphericalProduct_sampling(x_n, arclength);
        plot3(point(1, :), point(2, :), point(3, :), '.', 'Color', 'r')
        hold on
        plot3(point_fit(1, :), point_fit(2, :), point_fit(3, :), '.', 'Color', 'b')
        title(['Result after the ', num2str(iterEM), 'th EM iteration'])
        axis equal
        hold off
        disp(['Configuration after the ', num2str(iterEM), 'th EM iteration is:'])
        disp(x_n)
        disp(['Previous sigma2 = ', num2str(sigma2), ', Updated sigma2 = ', num2str(sigma2_n)])
        disp(['Previous Cost = ', num2str(cost), ', Updated Cost = ', num2str(cost_n), ', relative cost change = ', num2str(relative_cost)])
        disp('---------------------------------------------------------------')
        waitforbuttonpress
    end
    %--------------------------------------------------------------------------
    
    if (cost_n < toleranceEM && iterEM > 1) || (relative_cost < relative_toleranceEM && switched >= para.max_switch && iterEM > iterEM_min) % >1 >5
        cost = cost_n;
        sigma2 = sigma2_n;
        x = x_n;
        break
    end
    if relative_cost < relative_toleranceEM && iterEM ~= 1 % set different tolerance for switch and termination
        % activate switching algorithm to avoid local minimum
        switch_success = 0;
        % case1 - axis-mismatch similarity
        axis_0 = eul2rotm(x(6 : 8));
        axis_1 = circshift(axis_0, [0, 2]);
        axis_2 = circshift(axis_0, [0, 1]);
        eul_1 = rotm2eul(axis_1);
        eul_2 = rotm2eul(axis_2);
        x_axis = [x(2), x(1), x(4), x(5), x(3), eul_1, x(9 : 11); ...
            x(2), x(1), x(5), x(3), x(4), eul_2, x(9 : 11)];
        % case2 - duality similarity
        scale_ratio = circshift(x(3 : 5), 2) ./ x(3 : 5);
        scale_idx = find(and(scale_ratio > 0.9, scale_ratio < 1.1)); %0.5 1.5
        x_rot = zeros(size(scale_idx, 2), 11);
        rot_idx = 1;
        if ismember(1, scale_idx)
            eul_rot = rotm2eul(axis_0 * rotz(45));
            if x(2) <= 1
                x_rot(rot_idx, :) = [x(1), 2 - x(2), ((1 - sqrt(2)) * x(2) + sqrt(2)) * min(x(3), x(4)) * ones(1, 2), x(5), eul_rot, x(9 : 11)];
            else
                x_rot(rot_idx, :) = [x(1), 2 - x(2), ((sqrt(2)/2 - 1) * x(2) + 2 - sqrt(2)/2) * min(x(3), x(4)) * ones(1, 2), x(5), eul_rot, x(9 : 11)];
            end            
            rot_idx = rot_idx + 1;
        end
        if ismember(2, scale_idx)
            eul_rot = rotm2eul(axis_1 * rotz(45));
            if x(1) <= 1
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(4), x(5)) * ones(1, 2), x(3), eul_rot, x(9 : 11)];
            else
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((sqrt(2)/2 - 1) * x(1) + 2 - sqrt(2)/2) * min(x(4), x(5)) * ones(1, 2), x(3), eul_rot, x(9 : 11)];
            end    
            rot_idx = rot_idx + 1;
        end
        if ismember(3, scale_idx)
            eul_rot = rotm2eul(axis_2 * rotz(45));       
            if x(1) <= 1
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(5), x(3)) * ones(1, 2), x(4), eul_rot, x(9 : 11)];
            else
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((sqrt(2)/2 - 1) * x(1) + 2 - sqrt(2)/2) * min(x(5), x(3)) * ones(1, 2), x(4), eul_rot, x(9 : 11)];
            end    
        end
        
        % generate candidate configuration list with cost
        x_candidate = [x_axis; x_rot];
        cost_candidate = weighted_cost_switch(x_candidate, point, p);

        idx_nan = find(and(~isnan(cost_candidate), ~isinf(cost_candidate)));
        cost_candidate = cost_candidate(idx_nan);
        x_candidate = x_candidate(idx_nan, :);
        
        %------------------------------debug block---------------------------------
        if debug_mode == 1
            disp('ENTERING SWITCHING PROCESS')
            for i_debug = 1 : size(x_candidate, 1)
                figure
                [point_fit] = sphericalProduct_sampling(x_candidate(i_debug, :), arclength);
                plot3(point(1, :), point(2, :), point(3, :), '.', 'Color', 'r')
                hold on
                plot3(point_fit(1, :), point_fit(2, :), point_fit(3, :), '.', 'Color', 'b')
                title(['Candidate configuration ', num2str(i_debug)])
                axis equal
                hold off
                disp(['Candidate Configuration ', num2str(i_debug), ' is:'])
                disp(x_candidate(i_debug, :))
                disp(['Candidate Configuration ', num2str(i_debug), ' weighted cost is:'])
                disp(cost_candidate(i_debug))
                disp('---------------------------------------------------------------')
                waitforbuttonpress
            end
        end
        %--------------------------------------------------------------------------
        [~, idx] = sort(cost_candidate);
        for i_candidate = 1 : size(idx, 1)
                % adaptive upper bound
            if  adaptive_upper == 1
                R_current = eul2rotm(x_candidate(idx(i_candidate), 6 : 8));
                point_rot_current = R_current' * point - R_current' * x_candidate(idx(i_candidate), 9 : 11)';
%                 point_rot_current = R_current' * point;
                ub_a = 1.1 * [max(abs(point_rot_current(1, :))), max(abs(point_rot_current(2, :))), max(abs(point_rot_current(3, :)))];
                ub = [2.0 2.0 ub_a  2*pi 2*pi 2*pi ub_a];
                lb = [0.0 0.0 0.01 0.01 0.01 -2*pi -2*pi -2*pi -ub_a];
            end
            [x_switch, cost_switch] = lsqnonlin(cost_func, x_candidate(idx(i_candidate), :), lb, ub, options);
            if cost_switch < min(cost_n, cost)
                x = x_switch;
                cost = cost_switch;
                % update sigma
                sigma2 = cost_switch / (3 * sum(p));
                switch_success = 1;
                %------------------------------debug block---------------------------------
                if debug_mode == 1
                    disp('SWITCHED SUCCESS')
                    figure
                    [point_fit] = sphericalProduct_sampling(x, arclength);
                    plot3(point(1, :), point(2, :), point(3, :), '.', 'Color', 'r')
                    hold on
                    plot3(point_fit(1, :), point_fit(2, :), point_fit(3, :), '.', 'Color', 'b')
                    title('Switched Result')
                    axis equal
                    hold off
                    disp('Switched Configuration is:')
                    disp(x)
                    disp('Switched Configuration weighted cost is:')
                    disp(cost)
                    disp('---------------------------------------------------------------')
                    waitforbuttonpress
                end
                %--------------------------------------------------------------------------
                break
            end
        end
        if switch_success == 0
            cost = cost_n;
            sigma2 = sigma2_n;
            x = x_n;
        end
        switched = switched + 1;
    else
        cost = cost_n;
        sigma2 = sigma2_n;
        x = x_n;
    end
end

% scale back
if para.rescale == 1
    x(3 : 5) = x(3 : 5) * scale;
    x(9 : 11) = x(9 : 11) * scale;
end

% transform back from the center of mass
x(9 : 11) = x(9 : 11) + t0';

%% Functions
% ------------------eigen analysis-----------------------------------------
    function [EigenVector, EigenValue] = EigenAnalysis(point)
        CovM = point * point' ./ size(point, 2);
        [EigenVector, EigenValue] = eig(CovM);
        EigenVector = flip(EigenVector, 2);
    end

% ------------------distance function--------------------------------------
    function [dist] = distance(X, para)
        % transform pose parameters into R matrix and t vector
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11);
        % align the point cloud to the superquadrics coordinate
        X_c = R' * X - R' * t';
        % calulate the radial distance of each point
        r_0 = vecnorm(X_c);
        dist = r_0 .* abs(((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
            ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ (-para(1) / 2) - 1);
    end

% ------------------correspondence calculation ----------------------------
    function [p] = correspendence(dist, sigma2, w, p0)
        c = (2 * pi * sigma2) ^ (- 3 / 2);
        const = (w * p0) / (c * (1 - w));
        p = exp(-1 / (2 * sigma2) * dist .^ 2);
        p = p ./ (const + p);
    end

% ------------------weighed distance function -----------------------------
    function [value] = weighted_dist(para, X, p, sigma2)
        if sigma2 > 1e-10
            value = p .^ (1 / 2) .* distance(X, para);
        else
            value = abs((p .* distance(X, para) .^ 2 + 2 * sigma2 * log(surface_area(para)))) .^ (1 / 2); %0.0002
        end

    end

% ------------------surface area function ---------------------------------
    function [area] = surface_area(para)
        a00 = 8 * (para(3)*para(4) + para(4)*para(5) + para(3)*para(5)); % surface area of the cuboid
        a02 = 8 * (para(3)^2 + para(4)^2)^(1/2) * para(5) + 2 * para(3) * para(4);
        a20 = 4 * (para(3)*(para(4)^2 + para(5)^2)^(1/2) + para(4)*(para(3)^2 + para(5)^2)^(1/2));
        a = (para(3)^2 + para(4)^2)^(1/2);
        b = (para(4)^2 + para(5)^2)^(1/2);
        c = (para(3)^2 + para(5)^2)^(1/2);
        s = (a + b + c) / 2;
        a22 = 8 * (s * (s - a) * (s - b) * (s -c)) ^ (1/2);
        area = [1 - para(1) / 2, para(1) / 2] * [a00, a02; a20, a22] * [1 - para(2) / 2; para(2) / 2];
    end   
        
% ------------------weighed distance function (switch)---------------------
    function [value] = weighted_cost_switch(para, X, p)
        value = zeros(size(para, 1), 1);
        for i = 1 : size(para, 1)
            value(i, 1) = sum(p .* (distance(X, para(i, :)) .^ 2), 2);
        end
    end

end