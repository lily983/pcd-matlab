function results = getResultsR3(sampleNumber, objectType, Sigmax, sizeScale)
% getResultsR3 Generate random objects s1 and s2 and return PCD value
% generate by different methods
%
% Inputs
%     sizeScale                 object size scale
%     sampleNumber      number of random samples
%     objectType              object is sphere or ellipsoid
%
% Outputs
%     results                           struct array variables, including
%     collision status and PCD results by each methos

% Set default sizeScale to 0.1, object's size range from [0.1, 1)
if nargin == 3
    sizeScale = 0.1;
end

for iSample=1:sampleNumber
    %     Generate objects with random size, position, quaternion
    if strcmp(objectType, 'sphere')
        s1 = SuperQuadrics({(rand(1)+0.1)*10*sizeScale*ones(1,3), [1,1], [0, 0]...
            rand(3,1)*sizeScale, getRandomQuaternion(), [10, 10]});
        s2 = SuperQuadrics({(rand(1)+0.1)*10*sizeScale*ones(1,3),  [1,1], [0, 0]...
            (rand(3,1)+0.3)*10*sizeScale, getRandomQuaternion(),[10, 10]});
    elseif strcmp(objectType, 'ellipsoid')
        s1 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale, [1,1], [0, 0]...
            rand(3,1)*sizeScale, getRandomQuaternion(), [10, 10]});
        s2 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale,  [1,1], [0, 0]...
            (rand(3,1)+0.3)*10*sizeScale, getRandomQuaternion(),[10, 10]});
    else
        error('Input object type must be sphere or ellipsoid!');
    end
    
    [flag, dist, ~, ~] = collision_cfc(s1,s2, 'constrained');
    
    flagArray(iSample) = round(flag);
    distArray(iSample) = dist;
    
    %     Get PCD values of each methods
    [PCDExactArray(iSample), exactT(iSample)] = getPCDR3(s1, s2, Sigmax, 'PCD-exact');
    [PCDConvexArray(iSample), convexT(iSample)]  =  getPCDR3(s1, s2, Sigmax, 'PCD-convex');
    [PCDEB99Array(iSample), EB99T(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-EB-99');
    [PCDEB95Array(iSample), EB95T(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-EB-95');
    [PCDMaxpdfArray(iSample), maxpdfT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-maxpdf');
    [PCDGMMArray(iSample), GMMT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-GMM');
    [PCDEllipExact(iSample), EllipExactT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-UB-exact');
    [PCDEllipApproximation(iSample), EllipApproximationT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-UB-approximation');
    
    s1Array(:,iSample) = [s1.a, s1.q, s1.tc'];
    s2Array(:,iSample) = [s2.a, s2.q, s2.tc'];
end

%sort data by PCD-exact in ascent way and store sorted data in a struct
%array format
[results.PCDExactArray, index] = sort(PCDExactArray(1,1:size(PCDExactArray,2)), 2);
results.flagArray = flagArray(index);
results.distArray = distArray(index);
results.PCDConvexArray = PCDConvexArray(index);
results.PCDEB99Array = PCDEB99Array(index);
results.PCDEB95Array = PCDEB95Array(index);
results.PCDMaxpdfArray = PCDMaxpdfArray(index);
results.PCDGMMArray = PCDGMMArray(index);
results.PCDEllipExactArray = PCDEllipExact(index);
results.PCDEllipApproximationArray = PCDEllipApproximation(index);

% Computation time
results.exactT = exactT;
results.convexT = convexT;
results.EB99T = EB99T;
results.EB95T = EB95T;
results.maxpdfT = maxpdfT;
results.GMMT = GMMT;
results.EllipExactT = EllipExactT;
results.ElliPApproximationT = EllipApproximationT;
% s1 s2 Information
results.s1Array = s1Array(:, index);
results.s2Array = s2Array(:, index);
end