function results = getResultsR3QuadratureFormPaper(sampleNumber, objectType)
% Replicate experiments setting in paper "Tight Collision Probability for
% UAV Motion Planning in Uncertain Environment"

for iSample=1:sampleNumber
    % Generate Sigmax (double position error)
    Sigmax = [(2-0.01)*rand(1)+0.01, 0, 0;...
        0, (2-0.01)*rand(1)+0.01, 0;...
        0, 0, (2-0.01)*rand(1)+0.01];
    
    %     Generate objects with random size, position, quaternion
    if strcmp(objectType, 'sphere')
        s1 = SuperQuadrics({ones(1,3).*((2-0.2)*rand(1)+0.2), [1,1], [0, 0]...
            [0;0;0], [1,0,0,0], [10, 10]});
        s2 = SuperQuadrics({ones(1,3).*((2-0.2)*rand(1)+0.2),  [1,1], [0, 0]...
            4*rand(3,1), getRandomQuaternion(),[10, 10]});
    elseif strcmp(objectType, 'ellipsoid')
        s1 = SuperQuadrics({(2-0.2).*rand(1,3)+0.2, [1,1], [0, 0]...
            [0;0;0], [1,0,0,0], [10, 10]});
        s2 = SuperQuadrics({(2-0.2).*rand(1,3)+0.2,  [1,1], [0, 0]...
            4*rand(3,1), getRandomQuaternion(),[10, 10]});
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
%     PCDEB99Array(iSample)=NaN;PCDEB95Array(iSample)=NaN;
    [PCDMaxpdfArray(iSample), maxpdfT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-maxpdf');
    [PCDGMM2SGArray(iSample), GMM2SGT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-GMM');
%     [PCDGMM4SGArray(iSample), GMM4SGT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-GMM-4SG');
%     [PCDGMM5SGArray(iSample), GMM5SGT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-GMM-5SG');
    [PCDEllipExact(iSample), EllipExactT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-UB-exact');
    [PCDEllipApproximation(iSample), EllipApproximationT(iSample)]  = getPCDR3(s1, s2, Sigmax, 'PCD-ellip-UB-approximation');
%     PCDEllipBound(iSample)=NaN;EllipBoundT(iSample)=0;
    
    s1Array(:,iSample) = [s1.a, s1.q, s1.tc'];
    s2Array(:,iSample) = [s2.a, s2.q, s2.tc'];
    
    SigmaxArray(:,iSample) = diag(Sigmax);
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
results.PCDGMM2SGArray = PCDGMM2SGArray(index);
% results.PCDGMM4SGArray = PCDGMM4SGArray(index);
% results.PCDGMM5SGArray = PCDGMM5SGArray(index);
results.PCDEllipExactArray = PCDEllipExact(index);
results.PCDEllipApproximationArray = PCDEllipApproximation(index);

% Computation time
results.exactT = exactT;
results.convexT = convexT;
% results.EB99T = EB99T;
% results.EB95T = EB95T;
results.maxpdfT = maxpdfT;
results.GMM2SGT = GMM2SGT;
% results.GMM5SGT = GMM5SGT;
results.EllipExactT = EllipExactT;
results.ElliPApproximationT = EllipApproximationT;
% s1 s2 Information
results.s1Array = s1Array(:, index);
results.s2Array = s2Array(:, index);

% Sigmax
results.Sigmax = SigmaxArray(:, index);
end