function test_gmm_3D(s1, s2)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

[a1,a2,b1,b2] = get_gmm_param(Sigmaf, 3, 3, 1.15);

[X, Y, Z] = meshgrid(-200:2:200);

n=1;
I = zeros(size(X,1),size(X,1),size(X,1));
M = zeros(size(X,1),size(X,1),size(X,1));

for i=1:1:size(X,1)
    for j=1:1:size(X,1)
        for k=1:1:size(X,1)
             I(i, j,k) =  indicator_ellip([X(i,j,k); Y(i,j,k); Z(i,j,k)], mf, Sigmaf);
             M(i,j,k) = (a1 * exp( -0.5 * b1 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) )) ...
                 + a2 * exp( -0.5 * b2 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) ))) ...
                 / ((2*pi)^1.5*sqrt(det(Sigmaf)));
             if M(i,j,k) <= I(i,j,k)
                 P(:,n) = X(i,j,k); Y(i,j,k); Z(i,j,k);
                 n=n+1
             end
        end
    end
end

if n~=1
    figure; 
    scatter3(P(1,:), P(2,:), P(3,:));
end

end
