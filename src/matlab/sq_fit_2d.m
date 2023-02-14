classdef sq_fit_2d < handle
    %SQ_FIT_2D Given N 2D points, fit SuperEllipse
    %
    % Author: Sipu Ruan, ruansp@jhu.edu, Johns Hopkins University, 2019
    
    properties
        a        % shape parameter
        eps      % exponential parameter
        pose     % pose info: pose.tc -- center,
        %            pose.th -- angle of rotation
        pt       % 2D points: N x 2 array
        pt_conv  % convex hull of points
        
        cost     % minimized function value
    end
    
    methods
        function obj = sq_fit_2d(pts)
            obj.pt = pts';
            k = convhull(obj.pt');
            obj.pt_conv = obj.pt(:,k);
        end
        
        %% Superquadric model
        function F = sq_2d_model(obj, lambda)
            % Objective function to be minimized
            
            % parameters
            a1 = lambda(1);
            a2 = lambda(2);
            eps1 = lambda(3);
            tc = lambda(4:5)';
            th = lambda(6);
            
            % transformations
            R = rot2(th);
            x = R' * (obj.pt_conv - tc);
            
            % objective
            f = abs(x(1,:)./a1).^(2/eps1) + abs(x(2,:)./a2).^(2/eps1);
            F = sum(sqrt(a1*a2)*(f.^(eps1)-1).^2);
        end
        
        function [a, tc, th] = estimate_size(obj)
            %             tc = mean(obj.pt,2);
            %
            %             [U,S,~] = eig((obj.pt-tc)*(obj.pt-tc)');
            %             a = diag(S).^(0.5)/20;
            %
            %             if det(U) < 0
            %                 U(:,2) = -U(:,2);
            %             end
            %             th = atan2(U(2,1), U(1,1));
            
            [A, tc] = MinVolEllipse(obj.pt_conv, 1e-4);
            [R, D] = eig(A);
            a = diag(D).^(-0.5);
            th = atan2(R(2,1), R(1,1));
            %
            %             tc = mean(obj.pt,2);
            %             th = 0;
            %             pts_centered = obj.pt-tc;
            %             a = max(sqrt(sum(pts_centered.^2,1))) .* ones(2,1);
        end
        
        function [c, ceq] = constraint(obj, lambda)
            % parameters
            a1 = lambda(1);
            a2 = lambda(2);
            eps1 = lambda(3);
            tc = lambda(4:5)';
            th = lambda(6);
            
            % transformations
            R = rot2(th);
            x = R' * (obj.pt_conv - tc);
            
            % objective
            f = abs(x(1,:)./a1).^(2/eps1) + abs(x(2,:)./a2).^(2/eps1);
            c = f-1;
            ceq = [];
        end
        
        %% Optimizations
        function optimize(obj, xinit)
            if isempty(xinit)
                x0 = rand(1,6);
            else
                x0 = xinit;
            end
            
            lb = [0,0,0.1,-inf,-inf,-pi];
            ub = [inf,inf,2,inf,inf,pi];
            nlcon = @(lambda) obj.constraint(lambda);
            
            [lambda, obj.cost] = fmincon(@(lambda) obj.sq_2d_model(lambda),...
                x0, [],[],[],[],lb,ub,nlcon);
            
            %             options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
            %             [lambda, obj.cost] = fsolve(@(lambda) obj.sq_2d_model(lambda),...
            %                 x0, options);
            
            % Retrieve results
            obj.a = lambda(1:2);
            obj.eps = lambda(3);
            obj.pose.tc = lambda(4:5)';
            obj.pose.th = lambda(6);
        end
    end
end