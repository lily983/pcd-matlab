function results = getResultsR3(sizeScale, sampleNumber, objectType)
% getResultsR3 Generate random objects s1 and s2 and return PCD value
% generate by different methods 
%
% Inputs
%     sizeScale                 object size scale. If sizeScale=1, object
%     size range from (1, 10) 
%     sampleNumber      number of random samples 
%     objectType              object is sphere or ellipsoid.
%
% Outputs
%     results                           struct array variables, including
%     collision status and PCD results by each methos

for iSample=1:sampleNumber
%     Generate objects with random size, position, quaternion 
    if strcmp(objectType, 'sphere')
        s1 = SuperQuadrics({(rand(1)+1)*1*sizeScale*ones(1,3), [1,1], [0, 0]...
            rand(3,1)*0.3*sizeScale, getRandomQuaternion(), [20, 20]});
        s2 = SuperQuadrics({(rand(1)+1)*1*sizeScale*ones(1,3),  [1,1], [0, 0]...
            rand(3,1)*3.6*sizeScale, getRandomQuaternion(),[20,20]});
    elseif strcmp(objectType, 'ellip')
        s1 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale, [1,1], [0, 0]...
            rand(3,1)*sizeScale, getRandomQuaternion(), [20, 20]});
        s2 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale,  [1,1], [0, 0]...
            (rand(3,1)+0.3)*10*sizeScale, getRandomQuaternion(),[20,20]});
    else
        error('Input object type must be sphere or ellipsoid!');
    end

    [flag, dist, ~, ~] = collision_cfc(s1,s2, 'constrained');
 
%     Covariance matrix of the position error distribution
    Sigmax = zeros(3);
    Sigmax(1,1) = 8.0000e-00;
    Sigmax(2,2) = 8.0000e-00;
    Sigmax(3,3) = 8.0000e-00;
    Sigmax = Sigmax.*0.01*sizeScale^2;

    flagArray(iSample) = round(flag);
    distArray(iSample) = dist;
    
%     Get PCD values of each methods
    PCDExactArray(iSample) = getPCDR3(s1, s2, Sigmax, 'PCD-exact');
    PCDConvexArray(iSample) =  getPCDR3(s1, s2, Sigmax, 'PCD-convex');
    PCDMaxpdfArray(iSample) = getPCDR3(s1, s2, Sigmax, 'PCD-maxpdf');
    PCDSGArray(iSample) = getPCDR3(s1, s2, Sigmax, 'PCD-SG');
    PCDGMMArray(iSample) = getPCDR3(s1, s2, Sigmax, 'PCD-GMM');
end

%sort data by PCD-exact in ascent way and store sorted data in a struct
%array format
[results.PCDExactArray, index] = sort(PCDExactArray(1,1:size(PCDExactArray,2)), 2);
results.flagArray = flagArray(index);
results.distArray = distArray(index);
results.PCDConvexArray = PCDConvexArray(index);
results.PCDMaxpdfArray = PCDMaxpdfArray(index);
results.PCDSGArray = PCDSGArray(index);
results.PCDGMMArray = PCDGMMArray(index);
end