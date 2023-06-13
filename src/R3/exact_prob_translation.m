function prob = exact_prob_translation(s1, s2, Sigma, N)
%Calculate the MSC approximated exact probability

%check dimension
dimension = size(Sigma, 1);

prob = 0;

samples = mvnrnd(zeros(3,1), Sigma, N);

if dimension == 2
    s3 = SuperEllipse([s2.a(1), s2.a(2), s2.eps, s2.taper...
        s2.tc(1), s2.tc(2), s2.ang, s2.N]);
elseif dimension == 3
    s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
        s2.tc, s2.q, s2.N});
end

for i=1:N
    s3.tc = s2.tc + samples(i, :)';
    
    if dimension == 2
        %if s1 and s2 are sphere
        if isequal(s1.eps, 1) && isequal(s2.eps, 1) && isequal(s1.a./s1.a(1), ones(1,2)) && isequal(s2.a./s2.a(1), ones(1,2))
            if norm(s1.tc - s3.tc) <= (s1.a(1)+s3.a(1))
                prob = prob+1;
            end
        else
            if collision_mesh(s1, s2)
                prob=prob+1;
            end
        end
    elseif dimension == 3
        %if s1 and s2 are sphere
        if isequal(s1.eps, ones(1, 2)) && isequal(s2.eps, ones(1, 2)) && isequal(s1.a./s1.a(1), ones(1,3)) && isequal(s2.a./s2.a(1), ones(1,3))
            if norm(s1.tc - s3.tc) <= (s1.a(1)+s3.a(1))
                prob = prob+1;
            end
        else
            [flag, ~, ~, condition] = collision_cfc(s1,s3);
            if flag==1 && condition~=NaN
                prob = prob+1;
            elseif flag==1 && isnan(condition)
                [flag, ~, ~, condition] = collision_cfc(s1,s3, 'constrained');
                if isnan(condition)
                    prob=NaN;
                    break
                elseif flag==1
                    prob=prob+1;
                end
            end
        end
    end
    
end

prob = prob/N;

end