function   [prob, time, pdf_max_x, x_max] = maxPDFSQ(s1, s2, mx, Sigma)
% max_probability_pure_translation: Inspired by maxPDFSphere on
% superquadrics model (our method)
%
%Inputs
%   mx: Mean of relative position error x = x2-x1
%   Sigma: Covariance of the probability
%   s1, s2: SuperQuadrics
%Outputs
%   x_max: Center position of s2 that has the maximum probability of
%   colliding with s1
%   pdf_max_x: PDF at x_max
%   prob: Probability 
%   time: computation time
tic;
prob = 0;

% Function file shares the same pointer to s2 in the main file
s3 = SuperQuadrics({s2.a, s2.eps, [0, 0]...
    s2.tc, s2.q, s2.N});

if collision_cfc(s1, s2)
    x_max = s3.tc;
    pdf_max_x  = mvnpdf(x_max, mx, Sigma);
    prob = 1;
    time = toc;
    return
end

mink = MinkSumClosedForm(s1, s3, quat2rotm(s1.q), quat2rotm(s3.q));

R2 = euclideanfactory(2);

problem.M = R2;
problem.cost = @(phi)cost_function(phi, s1, s3, mink, Sigma);

% Solve the problem
options.maxiter = 100;
options.tolgradnorm = 1e-2;

[phi_opt, f_opt, ~] = trustregions(problem, [], options);

m1_opt = s1.GetGradientsFromSpherical(phi_opt');
x_max = mink.GetMinkSumFromGradient(m1_opt)  + s1.tc;

pdf_max_x = 1/( 8*pi^3 * sqrt(det(Sigma)) ) * exp(-1/2 * f_opt);

m1 = s1.GetGradientsCanonical();
x_mink = mink.GetMinkSumFromGradient(m1)+s1.tc;
SN = s1.N;
% X_ = reshape(x_mink(1, :), SN); 
% Y_= reshape(x_mink(2, :),  SN); 
% Z_ = reshape(x_mink(3, :),  SN); 
% patch_mink = surf2patch(X_,Y_,Z_, 'triangles');

alphaShape_mink = alphaShape(x_mink(1, :)', x_mink(2, :)', x_mink(3, :)', 3);
MinkowskiSumVolume = volume(alphaShape_mink);

prob = pdf_max_x * MinkowskiSumVolume;
time = toc;
end

function f = cost_function(phi, s1, s3, mink, Sigma)

m1 = s1.GetGradientsFromSpherical(phi');
mink_point = mink.GetMinkSumFromGradient(m1);

f = (mink_point - s3.tc)' / Sigma * (mink_point - s3.tc);

end