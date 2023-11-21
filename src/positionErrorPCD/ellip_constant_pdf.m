function prob=ellip_constant_pdf(s1, s2, mx, Sigmax, methodOptions)
% First do coordnate transformation to make the bounding ellipsoid as a
% unit ball located at the origin 
[mf, Sigmaf] = get_bounding_ellip(s1, s2);

new_mx = (sqrtm(Sigmaf)) \ (mx - mf);
new_Sigmax =  (sqrtm(Sigmaf)) \  Sigmax / sqrtm(Sigmaf);

if collision_cfc(s1,s2,'least-squares')
    prob=1;
    return
end

switch methodOptions
%     case 'pointA'
%         pdf = mvnpdf(zeros(1,3), new_mx', new_Sigmax);
    case 'pointB'
        pdf = mvnpdf(new_mx'./norm(new_mx), new_mx', new_Sigmax);
    case 'pointA'
        problem.M = spherefactory(3);
        problem.cost = @(x)(x-new_mx)'/new_Sigmax*(x-new_mx);
        options.maxiter = 100;
        options.tolgradnorm = 1e-2;
        x0=new_mx./norm(new_mx);
        [xopt, cost] = trustregions(problem, x0, options);
        pdf = mvnpdf(xopt'./norm(xopt), new_mx', new_Sigmax);
end

prob = pdf * pi * 4/3;

end