function prob = getPCDR3(s1, s2, Sigmax, method)
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
%           'PCD-SG': Upper bound of PCD-exact using single Gaussian
%           approximation for ellipsoids
%           'PCD-GMM': Upper bound of PCD-exact using Gaussian mixture
%           functuion for ellipsoids
%   Outputs:
%       prob probability of collision

%Check if s1 and s2 are sphere or ellipsoid
if isequal(s1.eps, ones(1, 2)) && isequal(s2.eps, ones(1, 2)) && isequal(s1.a./s1.a(1), ones(1,3)) && isequal(s2.a./s2.a(1), ones(1,3))
    objectType='sphere';
elseif isequal(s1.eps, ones(1, 2)) && isequal(s2.eps, ones(1, 2))
    objectType='ellipsoid';
end

switch method
    case 'PCD-exact'
        prob = exact_prob_translation(s1, s2, Sigmax, 1e+03);
    case 'PCD-maxpdf'
        if STRCMP(objectType, 'sphere')==1
            warning('Input objects are not sphere, unable to use PCD-maxpdf');
            prob = NaN;
            return
        end
        sphereVolume =  4*pi/3*(s1.a(1)+s2.a(1))^3 ;
        maxpdf = max_contact_probability_pure_translation(s2.tc, Sigmax, s1, s2);
        prob = sphereVolume*maxpdf;
    case 'PCD-convex'
        prob = max_contact_probability_polyhedron_pure_translation(Sigmax, s1, s2);
    case 'PCD-SG'
        if STRCMP(objectType, 'ellipsoid')==1
            warning('Input objects are not ellipsoid, unable to use PCD-SG');
            prob = NaN;
            return
        end
        prob = max_prob_single_gaussian(s1, s2, zeros(3,1), Sigmax);
    case 'PCD-GMM'
        if STRCMP(objectType, 'ellipsoid')==1
            warning('Input objects are not ellipsoid, unable to use PCD-GMM');
            prob = NaN;
            return
        end
        prob = max_prob_gaussian_mixture(s1, s2, mu, Sigma);
end

end