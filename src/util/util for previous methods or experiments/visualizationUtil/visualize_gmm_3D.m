function visualize_gmm_3D(s1, s2)
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

% [a1,a2,b1,b2] = get_gmm_param(Sigmaf, 3, 3, 1.15);
[c1,c2,c3,d1,d2,d3] = get_gmm_param_three(Sigmaf);
[a1,a2,a3, a4,b1,b2,b3,b4]=get_gmm_param_four(Sigmaf);

[X, Y, Z] = meshgrid(-5:0.1:5);

n=1;
I = zeros(size(X,1),size(X,1),size(X,1));
M1 = zeros(size(X,1),size(X,1),size(X,1));
M2 = zeros(size(X,1),size(X,1),size(X,1));

for i=1:1:size(X,1)
    for j=1:1:size(X,1)
        for k=1:1:size(X,1)
%              I(i, j,k) =  indicator_ellip([X(i,j,k); Y(i,j,k); Z(i,j,k)], mf, Sigmaf);
             M1(i,j,k) = (a1 * exp( -0.5 * b1 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) )) ...
                 + a2 * exp( -0.5 * b2 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) ))...
                 + a3 * exp( -0.5 * b3 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) ))...
                 + a4 * exp( -0.5 * b4 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) ))) ...
                 / ((2*pi)^1.5*sqrt(det(Sigmaf)));

            M2(i,j,k) = (c1 * exp( -0.5 * d1 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) )) ...
                 + c2 * exp( -0.5 * d2 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) ))...
                 + c3 * exp( -0.5 * d3 * ( ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf)' / Sigmaf * ([X(i,j,k); Y(i,j,k); Z(i,j,k)]-mf) ))) ...
                 / ((2*pi)^1.5*sqrt(det(Sigmaf)));
%              if M(i,j,k) <= I(i,j,k)
%                  n=n+1;
%                  result()M(i,j,k)-I(i,j,k)
%              end
%              result(n)=M1(i,j,k)-I(i,j,k);
             difference(n) = M1(i,j,k) -M2(i,j,k) ;
             n=n+1;
        end
    end
end
figure;
% plot(result)
plot(difference)

end
