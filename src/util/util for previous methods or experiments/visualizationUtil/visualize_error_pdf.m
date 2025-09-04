function visualize_error_pdf(mu, Sigma)
[X, Y] = meshgrid(-20:0.1:20);


for i=1:1:size(X,1)
    for j=1:1:size(X,1)
        PDF(i,j) = 1 * exp( -0.5 * ( ([X(i,j); Y(i,j)]-mu)' / Sigma * ([X(i,j); Y(i,j)]-mu) ))  / (2*pi*sqrt(det(Sigma)));
    end
end

sm = surf(X, Y, PDF);
sm.EdgeColor = 'none';
sm.FaceColor = 'interp';
sm.FaceAlpha = 0.8;
end