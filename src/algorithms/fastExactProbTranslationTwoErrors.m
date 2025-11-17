function [prob, time] = fastExactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, sampleNumbers)
%fastExactProbTranslationTwoErrors: this function
%is an implementation of paper 'A fast Monte Carlo algorithm for collision
%probability estimation'
%
%Inputs
%   s1, s2: Two convex objects (spheres, ellipsoids, superquadrics)
%   Sigma1: Covariance of the position error of object s1
%   Sigma2: Covariance of the position error of object s2
%   sampleNumbers: The number of Monte-Carlo sampling points
%Outputs
%   prob: The probability approximation
%   time: computation time

dimension = size(Sigma1, 1);

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

tic;
prob = 0;

% Here we sample the center of s1 and s2 based on the distribution of their
% position errors
samples1 = mvnrnd(s1.tc', Sigma1, sampleNumbers);
samples2 = mvnrnd(s2.tc', Sigma2, sampleNumbers);

% Here we create new objects mimicking s1 and s2, it is to prevent the s1
% s2 in the input get polluted inside the function because matlab gives the
% object handle
if dimension == 2
    s4 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
    s3 = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
        s1.tc(1), s1.tc(2), s1.ang, s1.N]);
elseif dimension == 3
    s4 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
        s2.tc, s2.q, s2.N});
    s3 = SuperQuadrics({s1.a, s1.eps, [0, 0]...
        s1.tc, s1.q, s1.N});
end

if dimension == 2
    %if s1 and s2 are sphere
    if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
        for i=1:sampleNumbers
            s3.tc = samples1(i, :)';
            s4.tc = samples2(i, :)';
            if norm(s3.tc - s4.tc) <= (s3.a(1)+s4.a(1))
                prob = prob+1;
            end
        end
    else
        for i=1:sampleNumbers
            s3.tc = samples1(i, :)';
            s4.tc = samples2(i, :)';
            if collision_GJK(s3, s4)
                prob=prob+1;
            end
        end
    end
elseif dimension == 3
    %if s1 and s2 are sphere
    if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
        for i=1:sampleNumbers
            s3.tc = samples1(i, :)';
            s4.tc = samples2(i, :)';
            if norm(s3.tc - s4.tc) <= (s3.a(1)+s4.a(1))
                prob = prob+1;
            end
        end
    elseif strcmp(s1objectType, 'ellip')==1 && strcmp(s2objectType, 'ellip')==1
        for i=1:sampleNumbers
            s3.tc = samples1(i, :)';
            s4.tc = samples2(i, :)';
            if collision_ellipsoid_asc(s3, s4)
                prob = prob+1;
            end
        end
     elseif strcmp(s1objectType, 'superquadrics')==1 && strcmp(s2objectType, 'superquadrics')==1
         for i=1:sampleNumbers
            s3.tc = samples1(i, :)';
            s4.tc = samples2(i, :)';
            % Here we do a double checking, because we use CFC and the
            % default solver uses fix-point-iteration, which computes
            % faster but often justify two seperated objects as collide, so
            % we use least-squares to double check if the two objects
            % collide or not
            if collision_cfc(s1, s4)
                if collision_cfc(s1, s4, 'least-squares')
                     prob = prob+1;
                end
            end
         end
    end
end


prob = prob/sampleNumbers;
time = toc;

end