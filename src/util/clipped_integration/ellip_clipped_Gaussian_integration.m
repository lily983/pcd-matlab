function r = ellip_clipped_Gaussian_integration(a, dimension)
if dimension==2
    r=first_error_integral(a);
elseif dimension==3
    r=sqrt(2/pi)*first_error_integral(a);
end
end