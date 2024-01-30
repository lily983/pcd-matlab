function [prob, t] = exactProbTranslation(s1, s2, Sigma, N)
%exact_prob_translation Calculate the MSC approximated exact probability
% Only for s2 with position error cases!
% This function is not fast MC for two error cases!
tic;
prob = 0;

%check dimension
dimension = size(Sigma, 1);

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

samples = mvnrnd(s2.tc, Sigma, N);

if dimension == 2
    s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
elseif dimension == 3
    s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
        s2.tc, s2.q, s2.N});
end

for i=1:N
    s3.tc = samples(i, :)';
    
    if dimension == 2
        %if s1 and s2 are sphere
        if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
            if norm(s1.tc - s3.tc) <= (s1.a(1)+s3.a(1))
                prob = prob+1;
            end
        else
            if collision_GJK(s1, s2)
                prob=prob+1;
            end
        end
    elseif dimension == 3
        %if s1 and s2 are sphere
        if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
            if norm(s1.tc - s3.tc) <= (s1.a(1)+s3.a(1))
                prob = prob+1;
            end
        elseif strcmp(s1objectType, 'ellip')==1 && strcmp(s2objectType, 'ellip')==1
            if collision_ellipsoid_asc(s1, s3)
                prob = prob+1;
            end
        elseif strcmp(s1objectType, 'superquadrics')==1 && strcmp(s2objectType, 'superquadrics')==1
            % Fix point iteration
            [flag, ~, ~, condition]=collision_cfc(s1, s3);
            % if condition doesn't satisfied, use constrained optimization
            % to compute again
            if isnan(condition) || condition>1e-03
                [flag, ~, ~, condition]=collision_cfc(s1, s3,'constrained');
                if ~isnan(condition) && flag
                    prob = prob+1;
%                     figure; hold on
%                     s1.PlotShape('b', 0.7);
%                     s2.PlotShape('g', 0.7);
%                     s3.PlotShape('m', 0.7);
                end
            end
        end
    end

end

prob = prob/N;
t = toc;
end