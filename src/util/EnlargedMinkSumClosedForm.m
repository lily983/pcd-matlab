classdef EnlargedMinkSumClosedForm < handle
    % EnlargedMinkSumClosedForm Computes exact closed-form Minkowski sums between 
    % two EnlargedSuperQuadrics
    %
    %  Inputs:
    %    sub1, sub2: EnlargedSuperQuadrics
    % 
    
    %% Variables
    properties
        sub1    % Class of convex smooth body
        sub2    % Class of convex smooth body
    end
 
    
    %% Functions
    methods
        %% Constructor
        function obj = EnlargedMinkSumClosedForm(sub1, sub2)
            obj.sub1 = sub1;
            obj.sub2 = sub2;
        end
        
        %% GetMinkSumFromNormal compute Minkowski sums from normals
        % Inputs:
        %   n  : normal vectors on EnlargedMinkSumClosedForm
        %
        % Output:
        %   mink: boundary points of EnlargedMinkSumClosedForm
        function mink = GetMinkSumFromNormal(obj, n)
            f_n1 = obj.sub1.GetPointsFromNormal( n );
            f_n2 = obj.sub2.GetPointsFromNormal( -n );
            
            mink = f_n1 - f_n2;
        end
        
        function mink = GetPoints(obj)
            n_1 = obj.sub1.sq.GetNormals();
            n_UB_1 = obj.sub1.R(:,:,1) * n_1;

            mink = obj.GetMinkSumFromNormal(n_UB_1);
        end
        
        function PlotShape(obj, color, faceAlpha, edgeAlpha)
            n_1 = obj.sub1.sq.GetNormals();
            n_UB_1 = obj.sub1.R(:,:,1) * n_1;

            x_mink = obj.GetMinkSumFromNormal(n_UB_1);
            N = obj.sub1.sq.N;
            x = reshape(x_mink(1,:), N(1), N(2));
            y = reshape(x_mink(2,:), N(1), N(2));
            z = reshape(x_mink(3,:), N(1), N(2));
            
            surf(x, y, z, 'FaceColor', color, 'EdgeColor', color,...
                'FaceAlpha', faceAlpha, 'EdgeAlpha', edgeAlpha);
        end
    end
end
