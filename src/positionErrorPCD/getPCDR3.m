function [prob, time] = getPCDR3(s1, s2, Sigmax, method)
%getPCDR3 Get the probability of collision of object s1 and s2 where s2
%subjects to the position error which follows the Gaussian distribution N(x; 0, Sigmax)
%   Inputs:
%       s1  object with exact pose
%       s2  object subjects to position error
%       Sigmax  covarian of position error distribution
%       method  name of method.
%           'PCD-exact': Monte-cale sampling approximation of exact
%           probability of collision
%           'PCD-maxpdf': Upper bound of PCD-exact using maximum pdf value
%           for spherical objects
%           'PCD-convex': Upper bound of PCD-exact using divergence theorm
%           for convex objects
%           'PCD-EB-99': Upper bound of  PCD-exact using enlarged bounding
%           object to approximate s2 (confidence level 99%)
%           'PCD-EB-95': Upper bound of  PCD-exact using enlarged bounding
%           object to approximate s2 (confidence level 95%)
%           'PCD-ellip-exact': Exact PCD of ellipsoids 
%           'PCD-ellip-bound': Upper bound for PCD of ellipsoids 
%           'PCD-SG': Upper bound of PCD-exact using single Gaussian
%           approximation for ellipsoids
%           'PCD-GMM': Upper bound of PCD-exact using Gaussian mixture
%           functuion for ellipsoids
%   Outputs:
%       prob probability of collision

switch method
    case 'PCD-exact'
        tic;
        prob = exact_prob_translation(s1, s2, Sigmax, 1e+03);
        time = toc;
    case 'PCD-maxpdf'
        tic;
        sphereVolume =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
        try 
            maxpdf = max_contact_probability_pure_translation(s2.tc, Sigmax, s1, s2);
            prob = sphereVolume*maxpdf;
        catch
            prob = NaN;
        end
        time = toc; 
    case 'PCD-convex'
        tic;
        try
            prob = max_contact_probability_polyhedron_pure_translation(Sigmax, s1, s2);
        catch 
            prob = NaN;
        end
        time = toc; 
    case 'PCD-EB-99'
        %This method is first write for SE3 case. For simplification, I
        %create the distribution for pose error which has exact rotation to
        %represent the distribution for position error
        tic;
        try
            muSE3 = [eye(3), s2.tc; 0, 0, 0, 1];
            SigmaSE3 = zeros(6,6); SigmaSE3(4:6,4:6) = Sigmax;
            prob = collision_enlarged_bounding_sq(s1, s2, muSE3, SigmaSE3,'99');
        catch 
            prob = NaN;
        end
        time = toc; 
        case 'PCD-EB-95'
        %This method is first write for SE3 case. For simplification, I
        %create the distribution for pose error which has exact rotation to
        %represent the distribution for position error
        tic;
        try
            muSE3 = [eye(3), s2.tc; 0, 0, 0, 1];
            SigmaSE3 = zeros(6,6); SigmaSE3(4:6,4:6) = Sigmax;
            prob = collision_enlarged_bounding_sq(s1, s2, muSE3, SigmaSE3, '95');
        catch 
            prob = NaN;
        end
        time = toc; 
    case 'PCD-ellip-exact'
        tic;
        try
            prob = max_and_exact_prob_ellipsoids(Sigmax, s1, s2, 'exact');
        catch 
            prob = NaN;
        end
        time = toc;
    case 'PCD-ellip-bound'
        tic;
        try
            prob = max_and_exact_prob_ellipsoids(Sigmax, s1, s2, 'max');
        catch 
            prob = NaN;
        end
        time = toc;
    case 'PCD-GMM-2SG'
        % Position error is assumed to follow N(x; 0, Sigmax)
        tic;
        try
            prob = max_prob_gaussian_mixture(s1, s2, zeros(3,1), Sigmax, 'accurate-result', '2SG');
        catch
            prob = NaN;
        end
        time = toc; 
    case 'PCD-GMM-4SG'
        tic;
        try
            prob = max_prob_gaussian_mixture(s1, s2, zeros(3,1), Sigmax, 'accurate-result', '4SG');
        catch
            prob = NaN;
        end
        time = toc; 
    case 'PCD-GMM-5SG'
        tic;
        try
            prob = max_prob_gaussian_mixture(s1, s2, zeros(3,1), Sigmax, 'accurate-result', '5SG');
        catch
            prob = NaN;
        end
        time = toc; 
end

% truncate the probability to 1 if it is larger than 1
if prob>1
    prob=1;
end

end