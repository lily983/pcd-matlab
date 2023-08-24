% function visualize_gmm_2D(mu, Sigma, bRatio, mb)
function visualize_gmm_2D(mu, Sigma)
% plot function value
[a1,a2,a3, a4,b1,b2,b3,b4]=get_gmm_param_four(Sigma);
[c1,c2,c3,d1,d2,d3] = get_gmm_param_three(Sigma);
[e1,e2,f1,f2] = get_gmm_param(Sigma, 10, 3);

% [X, Y] = meshgrid(round(mu)-15:0.05:round(mu)+15);
[X, Y] = meshgrid(-1:0.01:1);

Sigmax = eye(2, 2) * 8*1e-02;

for i=1:1:size(X,1)
    for j=1:1:size(X,1)
        M1(i,j) = a1*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b1)/b1 + ...
            a2*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b2)/b2 + ...
            a3*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b3)/b3 + ...
            a4*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b4)/b4;

        M2(i,j) = c1*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/d1)/d1 + ...
            c2*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/d2)/d2 + ...
            c3*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/d3)/d3;
        
        M3(i,j) = e1*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/f1)/f1 + ...
            e2*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/f2)/f2;

%         S(i,j) = mvnpdf([X(i,j), Y(i,j)], mu', Sigma) *  (2*pi)*sqrt(det(Sigma)) / exp(-1/2);

%         Sigmax = eye(2, 2) * 8*1e-02;
%         PDF(i,j) = mvnpdf([X(i,j), Y(i,j)], zeros(1,2), Sigmax);
        
        I(i, j) =  indicator_ellip([X(i,j); Y(i,j)], mu, Sigma);
        
    end
end

figure;
hold on;
% ss = surf(X, Y, S);
% ss.EdgeColor = 'none';
% ss.FaceColor = 'interp';
% ss.FaceAlpha = 0.2;

% figure;
% sm = surf(X, Y, M1);
% sm.EdgeColor = 'none';
% sm.FaceColor = 'interp';
% sm.FaceAlpha = 0.6;

sm2 = surf(X, Y, M2);
sm2.EdgeColor = 'none';
sm2.FaceColor = 'interp';
sm2.FaceAlpha = 0.6;

s3 = surf(X, Y, I);
s3.EdgeColor = 'none';
s3.FaceColor = 'interp';
s3.FaceAlpha = 0.6;

smn = surf(X, Y, M3);
smn.EdgeColor = 'none';
smn.FaceColor = 'interp';
smn.FaceAlpha = 0.7;

end


