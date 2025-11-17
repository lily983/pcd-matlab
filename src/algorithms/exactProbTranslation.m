function [prob, t] = exactProbTranslation(s1, s2, Sigma, sampleNumbers)
%exact_prob_translation: Monte-Carlo based numerical approximation of the
%collision probability between s1 and s2, where s2's position is subjected
%to Gaussian distributed errors
%
%Inputs
%   s1, s2: Two convex objects (spheres, ellipsoids, superquadrics)
%   Sigma: Covariance of the position error of object s2
%   sampleNumbers: The number of Monte-Carlo sampling points
%Outputs
%   prob: The probability approximation
%   time: computation time

%check dimension
dimension = size(Sigma, 1);

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

tic;
prob = 0;

% Here we sample the center of s2 based on the distribution of its
% position error
samples = mvnrnd(s2.tc, Sigma, sampleNumbers);

if dimension == 2
    s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
elseif dimension == 3
    s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
        s2.tc, s2.q, s2.N});
end

if dimension == 2
    %if s1 and s2 are sphere
    if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
        for i=1:sampleNumbers
            s3.tc = samples(i, :)';
            if norm(s1.tc - s3.tc) <= (s1.a(1)+s3.a(1))
                prob = prob+1;
            end
        end
    else
        for i=1:sampleNumbers
            s3.tc = samples(i, :)';
            if collision_GJK(s1, s2)
                prob=prob+1;
            end
        end
    end
elseif dimension == 3
    %if s1 and s2 are sphere
    if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
        for i=1:sampleNumbers
            s3.tc = samples(i, :)';
            if norm(s1.tc - s3.tc) <= (s1.a(1)+s3.a(1))
                prob = prob+1;
            end
        end
    elseif strcmp(s1objectType, 'ellip')==1 && strcmp(s2objectType, 'ellip')==1
        for i=1:sampleNumbers
            s3.tc = samples(i, :)';
            if collision_ellipsoid_asc(s1, s3)
                prob = prob+1;
            end
        end
    elseif strcmp(s1objectType, 'superquadrics')==1 && strcmp(s2objectType, 'superquadrics')==1
        for i=1:sampleNumbers
            s3.tc = samples(i, :)';
            % Fix point iteration
            [~, ~, ~, condition]=collision_cfc(s1, s3);
            % if condition doesn't satisfied, use constrained optimization
            % to compute again
            if isnan(condition) || condition>1e-03
                [flag, ~, ~, condition]=collision_cfc(s1, s3,'constrained');
                if ~isnan(condition) && flag
                    prob = prob+1;
                end
            end
        end
    end
end

end

prob = prob/N;
t = toc;
end