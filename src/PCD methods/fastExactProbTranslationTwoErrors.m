function [prob, time] = fastExactProbTranslationTwoErrors(s1, s2, Sigma1, Sigma2, sampleNumbers)
%fastExactProbTranslationTwoErrors: Fast Monte Carlo sampling-based
%approximation of the exact PCD of two position errors cases. This function
%is an implementation of paper 'A fast Monte Carlo algorithm for collision
%probability estimation'
tic;
prob = 0;

dimension = size(s1.a,2);

%Check if s1 and s2 are superquadric 
s1objectType = getObjectType(s1);
s2objectType = getObjectType(s2);

% samples1 = mvnrnd(s1.tc', Sigma1, sampleNumbers);
samples2 = mvnrnd(s2.tc', Sigma2, sampleNumbers);

if dimension == 2
    s4 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
%     s3 = SuperEllipse([s1.a(1), s1.a(2), s1.eps, s1.taper...
%         s1.tc(1), s1.tc(2), s1.ang, s1.N]);
elseif dimension == 3
    s4 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
        s2.tc, s2.q, s2.N});
%     s3 = SuperQuadrics({s1.a, s1.eps, [0, 0]...
%         s1.tc, s1.q, s1.N});
end

for i=1:sampleNumbers
%     s3.tc = samples1(i, :)';
   
    s4.tc = samples2(i, :)';

    if dimension == 2
        %if s1 and s2 are sphere
        if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
            if norm(s1.tc - s4.tc) <= (s1.a(1)+s4.a(1))
                prob = prob+1;
            end
        else
            if collision_GJK(s1, s4)
                prob=prob+1;
            end
        end
    elseif dimension == 3
        %if s1 and s2 are sphere
        if strcmp(s1objectType, 'sphere')==1 && strcmp(s2objectType, 'sphere')==1
            if norm(s1.tc - s4.tc) <= (s1.a(1)+s4.a(1))
                prob = prob+1;
            end
        elseif strcmp(s1objectType, 'ellip')==1 && strcmp(s2objectType, 'ellip')==1
            if collision_ellipsoid_asc(s1, s4)
                prob = prob+1;
            end
         elseif strcmp(s1objectType, 'superquadrics')==1 && strcmp(s2objectType, 'superquadrics')==1
            if collision_cfc(s1, s4)
                if collision_cfc(s1, s4, 'least-squares')
                     prob = prob+1;
%                     figure; hold on;axis equal;axis off
%                     s1.PlotShape(hex2rgb('9fc5e8'), 0.6);
%                     s4.PlotShape(hex2rgb('45AC59'), 0.8);
                end
            end
        end
    end
end

prob = prob/sampleNumbers;
time = toc;

end