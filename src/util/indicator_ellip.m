function flag = indicator_ellip(x, mu, Sigma)
% flag=1: point x is inside ellipsoid 
condition = ( (x-mu)'/Sigma*(x-mu) <=1 );
if condition
    flag=1;
else 
    flag=0;
end
end