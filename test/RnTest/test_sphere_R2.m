clear;clc;

N=500;
for i=1:N
    % Generate s1 and s2 with random aize and pose
    randSize = (rand(1)+1.2)*10*ones(1,2);
    randPosi = rand(1,2)*3;
    s1 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
    
    randSize = (rand(1)+1.2)*10*ones(1,2);
    randPosi = rand(1,2)*40;
    s2 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
   flag = collision_mesh(s1, s2);
    
%     figure; hold on;
%     s1.PlotShape('b', 0.4);
%     s2.PlotShape('g', 0.4);
%     pause(1)

    mu=zeros(2,1);
    
    Sigma = zeros(2);
    Sigma(1,1) = 1.0000e-00;
    Sigma(2,2) = 2.0000e-00;
%     Sigma   = Sigma.*0.1;
    
     try
         [prob_single_gaussian, weightedNorm] = max_prob_single_gaussian(s1, s2, mu, Sigma);
     catch 
         prob_single_gaussian = NaN;
     end
     
     try 
         prob_gaussian_mixture = max_prob_gaussian_mixture(s1, s2, mu, Sigma, 10, 10, 3);
     catch
         prob_gaussian_mixture = NaN;
     end
    
    try
        exact_prob_x = exact_prob_translation(s1, s2, mu, Sigma, 1e+03);
    catch
        exact_prob_x = NaN;
    end
        
    j=1; result(j, i) = double(flag);
    j=j+1; result(j, i) = exact_prob_x;
    j=j+1; result(j, i) = prob_single_gaussian;
    j=j+1; result(j, i) = prob_gaussian_mixture;

    j=j+2; result(j, i) = s1.a(1);
    j=j+1; result(j:j+1, i) = s1.tc;
    j=j+2; result(j, i) =s1.ang;
    j=j+2; result(j, i) = s2.a(1);
    j=j+1; result(j:j+1, i) = s2.tc;
    j=j+2; result(j, i) =s2.ang;
    
    j=j+2; result(j, i) = weightedNorm;

end    

%% size*0.1
clear;clc;

N=500;
for i=1:N
    % Generate s1 and s2 with random aize and pose
    randSize = (rand(1)+1.2)*ones(1,2) * 1 ;
    randPosi = rand(1,2)*0.5;
    s1 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
    
    randSize = (rand(1)+1.2)*ones(1,2) * 1 ;
    randPosi = rand(1,2)*5;
    s2 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
   flag = collision_mesh(s1, s2);
    
%     figure; hold on;
%     s1.PlotShape('b', 0.4);
%     s2.PlotShape('g', 0.4);
%     pause(1)

    mu=zeros(2,1);
    
    Sigma = zeros(2);
    Sigma(1,1) = 1.0000e-00;
    Sigma(2,2) = 2.0000e-00;
    Sigma = Sigma.*0.01;
    
     try
         [prob_single_gaussian, weightedNorm] = max_prob_single_gaussian(s1, s2, mu, Sigma);
     catch 
         prob_single_gaussian = NaN;
     end
     
     try 
         prob_gaussian_mixture = max_prob_gaussian_mixture(s1, s2, mu, Sigma, 10, 10, 3);
     catch
         prob_gaussian_mixture = NaN;
     end
    
    try
        exact_prob_x = exact_prob_translation(s1, s2, mu, Sigma, 1e+03);
    catch
        exact_prob_x = NaN;
    end
        
    j=1; result(j, i) = double(flag);
    j=j+1; result(j, i) = exact_prob_x;
    j=j+1; result(j, i) = prob_single_gaussian;
    j=j+1; result(j, i) = prob_gaussian_mixture;

    j=j+2; result(j, i) = s1.a(1);
    j=j+1; result(j:j+1, i) = s1.tc;
    j=j+2; result(j, i) =s1.ang;
    j=j+2; result(j, i) = s2.a(1);
    j=j+1; result(j:j+1, i) = s2.tc;
    j=j+2; result(j, i) =s2.ang;
    
    j=j+2; result(j, i) = weightedNorm;

end 

%% size*0.01
clear;clc;

N=500;
for i=1:N
    % Generate s1 and s2 with random aize and pose
    randSize = (rand(1)+1.2)*ones(1,2) * 0.1 ;
    randPosi = rand(1,2)*0.05;
    s1 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
    
    randSize = (rand(1)+1.2)*ones(1,2) * 0.1 ;
    randPosi = rand(1,2)*0.5;
    s2 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
   flag = collision_mesh(s1, s2);
    
%     figure; hold on;
%     s1.PlotShape('b', 0.4);
%     s2.PlotShape('g', 0.4);
%     pause(1)

    mu=zeros(2,1);
    
    Sigma = zeros(2);
    Sigma(1,1) = 1.0000e-00;
    Sigma(2,2) = 2.0000e-00;
    Sigma = Sigma.*0.0001;
    
     try
         [prob_single_gaussian, weightedNorm] = max_prob_single_gaussian(s1, s2, mu, Sigma);
     catch 
         prob_single_gaussian = NaN;
     end
     
     try 
         prob_gaussian_mixture = max_prob_gaussian_mixture(s1, s2, mu, Sigma, 10, 10, 3);
     catch
         prob_gaussian_mixture = NaN;
     end
    
    try
        exact_prob_x = exact_prob_translation(s1, s2, mu, Sigma, 1e+03);
    catch
        exact_prob_x = NaN;
    end
        
    j=1; result(j, i) = double(flag);
    j=j+1; result(j, i) = exact_prob_x;
    j=j+1; result(j, i) = prob_single_gaussian;
    j=j+1; result(j, i) = prob_gaussian_mixture;

    j=j+2; result(j, i) = s1.a(1);
    j=j+1; result(j:j+1, i) = s1.tc;
    j=j+2; result(j, i) =s1.ang;
    j=j+2; result(j, i) = s2.a(1);
    j=j+1; result(j:j+1, i) = s2.tc;
    j=j+2; result(j, i) =s2.ang;
    
    j=j+2; result(j, i) = weightedNorm;

end 

%% size*10
clear;clc;

N=500;
for i=1:N
    % Generate s1 and s2 with random aize and pose
    randSize = (rand(1)+1.2)*ones(1,2) * 100 ;
    randPosi = rand(1,2)*3;
    s1 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
    
    randSize = (rand(1)+1.2)*ones(1,2) * 100 ;
    randPosi = rand(1,2)*350;
    s2 = SuperEllipse([randSize(1), randSize(2), 1, 0 ...
        randPosi(1), randPosi(2), rand(1)*2*pi, 25]);
   flag = collision_mesh(s1, s2);
    
%     figure; hold on;
%     s1.PlotShape('b', 0.4);
%     s2.PlotShape('g', 0.4);
%     pause(1)

    mu=zeros(2,1);
    
    Sigma = zeros(2);
    Sigma(1,1) = 1.0000e-00;
    Sigma(2,2) = 2.0000e-00;
    Sigma = Sigma.*10;
    
     try
         [prob_single_gaussian, weightedNorm] = max_prob_single_gaussian(s1, s2, mu, Sigma);
     catch 
         prob_single_gaussian = NaN;
     end
     
     try 
         prob_gaussian_mixture = max_prob_gaussian_mixture(s1, s2, mu, Sigma, 10, 10, 3);
     catch
         prob_gaussian_mixture = NaN;
     end
    
    try
        exact_prob_x = exact_prob_translation(s1, s2, mu, Sigma, 1e+03);
    catch
        exact_prob_x = NaN;
    end
        
    j=1; result(j, i) = double(flag);
    j=j+1; result(j, i) = exact_prob_x;
    j=j+1; result(j, i) = prob_single_gaussian;
    j=j+1; result(j, i) = prob_gaussian_mixture;

    j=j+2; result(j, i) = s1.a(1);
    j=j+1; result(j:j+1, i) = s1.tc;
    j=j+2; result(j, i) =s1.ang;
    j=j+2; result(j, i) = s2.a(1);
    j=j+1; result(j:j+1, i) = s2.tc;
    j=j+2; result(j, i) =s2.ang;
    
    j=j+2; result(j, i) = weightedNorm;

end 

%% Plot results
exact = result(2, 1:N);
our_bound = result(3, 1:N);
our_bound_GMM = result(4, 1:N);
weighted_norm = result(16, 1:N);
flag = result(1, 1:N);

% [sort_exact, index]=sort(exact,2);
[sort_norm, index] = sort(weighted_norm, 'descend');

for i=1:N
    sort_flag(i) = flag(index(i));
    sort_exact(i) = exact(index(i));
    sort_our_bound(i) = our_bound(index(i));
    sort_our_bound_GMM(i) = our_bound_GMM(index(i));
end
%%
figure; hold on;
scatter(1:N, sort_exact(1:N), 'k','diamond');
scatter(1:N,  sort_our_bound(1:N),  'r', '+');
scatter(1:N,  sort_our_bound_GMM(1:N), 'b', '+');

%%
figure; hold on;
plot(1:N, sort_exact(1:N), 'k')
plot(1:N, sort_our_bound(1:N), 'r')
plot(1:N, sort_our_bound_GMM(1:N), 'b')
plot(1:N, sort_norm(1:N), 'm');
plot(1:N, sort_flag(1:N), 'c');
plot(1:N, ones(N, 1), 'y');

