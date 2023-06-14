function visualize_gmm_2D(mu, Sigma, aRatio, bRatio, mb)
% plot function value
[a1,a2,b1,b2] = get_gmm_param(Sigma, aRatio, bRatio, mb);

[X, Y] = meshgrid(round(mu)-100:1:round(mu)+100);


for i=1:1:size(X,1)
    for j=1:1:size(X,1)
        M(i,j) = a1*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b1)/b1 + a2*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b2)/b2;
%         M(i,j) = (a1 * exp( -0.5 * b1 * ( ([X(i,j); Y(i,j)]-mu)' / Sigma * ([X(i,j); Y(i,j)]-mu) )) ...
%             + a2 * exp( -0.5 * b2 * ( ([X(i,j); Y(i,j)]-mu)' / Sigma * ([X(i,j); Y(i,j)]-mu) ))) ...
%             / (2*pi*sqrt(det(Sigma)));
        
        S(i,j) = mvnpdf([X(i,j), Y(i,j)], mu', Sigma) *  (2*pi)*sqrt(det(Sigma)) / exp(-1/2);
        
        I(i, j) =  indicator_ellip([X(i,j); Y(i,j)], mu, Sigma);
        
%         DM(i, j) = M(i,j) - I(i,j);
%         DS(i,j) = S(i,j) - I(i, j);
        
%         N(i,j) = a1 * mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b1) / sqrt(b1) ...
%             + a2 * mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b2) / sqrt(b2);
%         DMN(i,j) = M(i,j) - N(i,j);
            
%             V=[X(i,j), Y(i,j)];
%             T1(i,j) = a1* exp(-0.5*b1*( (V'-mu)'/Sigma*(V'-mu) )) / (2*pi*sqrt(det(Sigma))) ;
% %             T2(i,j) = a2* exp(-0.5*(V'-mu)'/(Sigma/b2)*(V'-mu)) / (2*pi*sqrt(det(Sigma))) ;
%             T2(i,j) = -a2* exp(-0.5*b2*( (V'-mu)'/Sigma*(V'-mu) )) / (2*pi*sqrt(det(Sigma))) ;
%             DT(i,j) = -T1(i,j) + T2(i,j);
    end
end

% figure;
% surf(X, Y, Z);
% figure;
% surf(X, Y, G);


% figure;
% hold on;
ss = surf(X, Y, S);
ss.EdgeColor = 'none';
ss.FaceColor = 'interp';
ss.FaceAlpha = 0.3;

% % figure;
sm = surf(X, Y, M);
sm.EdgeColor = 'none';
sm.FaceColor = 'interp';
sm.FaceAlpha = 0.5;
% % 
% % 

s3 = surf(X, Y, I);
s3.EdgeColor = 'none';
s3.FaceColor = 'interp';
s3.FaceAlpha = 0.6;

% smn = surf(X, Y, DMN);
% smn.EdgeColor = 'none';
% smn.FaceColor = 'interp';
% smn.FaceAlpha = 0.3;

% sn = surf(X, Y, N);
% sn.EdgeColor = 'none';
% sn.FaceColor = 'interp';
% sn.FaceAlpha = 0.3;
% 
% figure
% st1 = surf(X, Y, T1);
% st1.EdgeColor = 'none';
% st1.FaceColor = 'interp';
% st1.FaceAlpha = 0.3;
% 
% % pause(2)
% st2 = surf(X, Y, T2);
% st2.EdgeColor = 'none';
% st2.FaceColor = 'interp';
% st2.FaceAlpha = 0.4;

% figure;
% sdt = surf(X, Y, DT);
% sdt.EdgeColor = 'none';
% sdt.FaceColor = 'interp';
% sdt.FaceAlpha = 0.3;
end

