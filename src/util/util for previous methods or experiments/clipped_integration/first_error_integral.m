function f1 = first_error_integral(a)
f1 = -a.*exp(-a.^2/2) + zero_error_integral(a);
end