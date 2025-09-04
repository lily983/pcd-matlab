function f0 = zero_error_integral(a)
%zero_error_integral Returns the value of zero_error_integral function at a
%defined for ellipsoid-clipped integration of a Gaussian distribution
%  Inputs:
%    a          The scale factor of integration region relative to a standard
%    ellipsoid defined by the mean and covariance of the Gaussian
%    distribution
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021

f0 = sqrt(pi/2)*erf(a/sqrt(2));
end