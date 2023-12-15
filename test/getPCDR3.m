function [prob, time] = getPCDR3(s1, s2, mx, Sigma1, Sigma2, method)
%getPCDR3 Get the probability of collision of object s1 and s2 where s2
%subjects to the position error which follows the Gaussian distribution N(x; 0, Sigmax)
%   Inputs:
%       s1  object with exact pose
%       s2  object subjects to position error
%       mx center of position error distribution
%       Sigmax  covarian of position error distribution
%       method  name of method.
%           'Exact': Monte-cale sampling approximation of exact
%           probability of collision (fast MC for two errors cases)
%           'Exact_two': Monte-cale sampling approximation of exact
%           probability of collision for two errors cases
%           probability of collision
%           'Maxpdf': Upper bound of PCD-exact using maximum pdf value
%           for spherical objects
%           'Maxpdf_SQ': (our method) Upper bound of PCD-exact using maximum pdf value
%           for superquadrics
%           'Divergence_mesh': Upper bound of PCD-exact using divergence theorm
%           for convex objects
%           'EB_99': Upper bound of  PCD-exact using enlarged bounding
%           object to approximate s2 (confidence level 99%)
%           'EB_95': Upper bound of  PCD-exact using enlarged bounding
%           object to approximate s2 (confidence level 95%)
%           'Quadratic_exact': Exact PCD of ellipsoids based on quadratic
%           form
%           'Quadratic_bound': Upper bound for PCD of ellipsoids based on
%           quadratic form
%           'LCC_center_point': Linear chance constraint, vector connecting
%           center of two objects
%           'LCC_center_point_cfc': Linear chance constraint, vector connecting
%           center of two objects ( using closed-form contact space without
%           the need to get bounding ellipsoid for Minkowski sum)
%            'LCC_closed_point': (our method): Linear chance constraint,
%            normal vector of the closed point on Minkowski sum  
%            'LCC_tangent_point': (our method) Linear chance constraint,
%            first transform position error Sigma to eye(3), then find
%            closed point on Minkowski sum to mx
%           'PCD-GMF': (our method) Upper bound of PCD-exact using Gaussian mixture
%           functuion for ellipsoids
%   Outputs:
%       prob probability of collision

Sigmax = Sigma1 + Sigma2;

switch method
    case 'Exact'
        [prob, time] = exactProbTranslation(s1, s2, Sigma2, 1e+03);
    case 'Exact_two'
        [prob, time] = exactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, 1e+03);
    case 'Fast_sampling'
        [prob, time] = fastExactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, 1e+03);
    case 'Maxpdf'
        try 
            [prob, time] = maxPDFSphere(s1, s2, mx, Sigma2);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'Maxpdf_SQ'
        try 
            [prob, time] = maxPDFSQ(s1, s2, mx, Sigma2);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'Divergence_mesh'
        try
            [prob, time] = divergenceMesh(s1, s2, mx, Sigma2);
        catch 
            prob = NaN;
            time = NaN;
        end
    case 'EB_99'
        try
            [prob, time] = enlargedBoundingVolume(s1, s2, Sigma1, Sigma2, 0.99);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'EB_95'
        try
            [prob, time] = enlargedBoundingVolume(s1, s2, Sigma1, Sigma2, 0.95);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'GMF'
        try 
            [prob, time] = gaussianMixturedFunction(s1, s2, mx, Sigmax);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'Quadratic_exact'
        try
            [prob, time] = quadraticExact(s1, s2, mx, Sigmax);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'Quadratic_bound'
        try
            [prob, time] = quadraticBound(s1, s2, mx, Sigmax);
        catch
            prob = NaN;
            time = NaN;
        end
    case 'LCC_tangent_point'
        try
            [prob, time] = linearChanceConstraintBound(s1, s2, mx, Sigmax, 'tangent-point', false);
        catch 
            prob = NaN;
            time = NaN;
        end
    case 'LCC_closed_point'
        try
            [prob, time] = linearChanceConstraintBound(s1, s2, mx, Sigmax, 'closed-point', false);
        catch 
            prob = NaN;
            time = NaN;
        end
    case 'LCC_center_point'
        try
            [prob, time] = linearChanceConstraintBound(s1, s2, mx, Sigmax, 'center-point', false);
        catch 
            prob = NaN;
            time = NaN;
        end
    case 'LCC_center_point_cfc'
        try
            [prob, time] = linearChanceConstraintBound(s1, s2, mx, Sigmax, 'center-point-cfc', false);
        catch 
            prob = NaN;
            time = NaN;
        end
end

% truncate the probability to 1 if it is larger than 1
if prob>1
    prob=1;
end

end