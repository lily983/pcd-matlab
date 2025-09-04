clc;clear; close all;

%%
% set mu and Sigma
mu = zeros(2,1);
Sigma =[1, 0; 0, 4];

s1 = SuperEllipse([1, 2, 1, 0, 0, 0, 0, 50]);

% plot function value
[a, b]=get_gmm_param(Sigma,'5SG');
[e,f] = get_gmm_param(Sigma, '2SG');

a1=a(1);a2=a(2);a3=a(3);a4=a(4);a5=a(5);
b1=b(1);b2=b(2);b3=b(3);b4=b(4);b5=b(5);

e1=e(1);e2 = e(2);
f1=f(1); f2=f(2);

[X, Y] = meshgrid(-2:0.01:2);

% Sigmax = eye(2, 2) * 8*1e-02;

for i=1:1:size(X,1)
    for j=1:1:size(X,1)
        M1(i,j) = a1*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b1)/b1 + ...
            a2*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b2)/b2 + ...
            a3*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b3)/b3 + ...
            a4*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b4)/b4 + ...
            a5*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/b5)/b5;
        
        M3(i,j) = e1*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/f1)/f1 + ...
            e2*mvnpdf([X(i,j), Y(i,j)], mu', Sigma/f2)/f2;

%         S(i,j) = mvnpdf([X(i,j), Y(i,j)], mu', Sigma) *  (2*pi)*sqrt(det(Sigma)) / exp(-1/2);

%         Sigmax = eye(2, 2) * 8*1e-02;
%         PDF(i,j) = mvnpdf([X(i,j), Y(i,j)], zeros(1,2), Sigmax);
        
        I(i, j) =  indicator_ellip([X(i,j); Y(i,j)], mu, Sigma);
        
    end
end

%%
figure; hold on;

sm2 = surf(X, Y, M1);
sm2.EdgeColor = 'none';
sm2.FaceColor = hex2rgb('F99B45');
sm2.FaceAlpha = 1;

s3 = surf(X, Y, I);
s3.EdgeColor = 'none';
s3.FaceColor = hex2rgb('D95980');
s3.FaceAlpha = 0.8;

smn = surf(X, Y, M3);
smn.EdgeColor = 'none';
smn.FaceColor = hex2rgb('284E60');
smn.FaceAlpha = 1;

zlim([0, 1])
ylim([0, 2])
xlim([0, 1.3])

s1_color = hex2rgb('D95980');
s1.PlotShape(s1_color, 0, 1);


