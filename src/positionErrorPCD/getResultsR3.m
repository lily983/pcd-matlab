function results = getResultsR3(sampleNumber, objectType, Sigmax, keySet)
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
if nargin == 4
    sizeScale = 0.1;
end

f = waitbar(0, 'Starting');

for iSample=1:sampleNumber
    waitbar(iSample/sampleNumber, f, sprintf('Progress: %d %%', floor(iSample/sampleNumber*100)));
    
    SN = [10, 10];
    
    if strcmp(objectType, 'sphere')
        s1 = SuperQuadrics({(rand(1)+0.1)*10*sizeScale*ones(1,3), [1,1], [0, 0]...
            rand(3,1)*sizeScale, getRandomQuaternion(), SN});
        s2 = SuperQuadrics({(rand(1)+0.1)*10*sizeScale*ones(1,3),  [1,1], [0, 0]...
            (rand(3,1)+0.3)*10*sizeScale, getRandomQuaternion(),SN});
    elseif strcmp(objectType, 'ellipsoid')
        s1 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale, [1,1], [0, 0]...
            rand(3,1)*sizeScale, getRandomQuaternion(), SN});
        s2 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale,  [1,1], [0, 0]...
            (rand(3,1)+0.3)*10*sizeScale, getRandomQuaternion(),SN});
    elseif strcmp(objectType, 'superquadrics')
        SN = [40, 40];
        s1 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale, 2*rand(1,2), [0, 0]...
            rand(3,1)*sizeScale, getRandomQuaternion(), SN});
        s2 = SuperQuadrics({(rand(1,3)+0.1)*10*sizeScale, 2*rand(1,2), [0, 0]...
            (rand(3,1)+0.3)*10*sizeScale, getRandomQuaternion(),SN});
    end
    
    [flag, dist, ~, condition] = collision_cfc(s1,s2, 'constrained');
    
    % Check cfc condition between s1 and s2, if not converge, then generate
    % new pairs
    if objectType=="superquadrics"
        if condition > 1e-02 || isnan(condition)
            iSample = iSample-1;
            continue;
        end
    end
    
    flagArray(iSample) = round(flag);
    distArray(iSample) = dist;
    
    mx = s2.tc-s1.tc;
    
    for i = 1:size(keySet,2)
        iMethod = string(keySet(i));
        [results.(genvarname(iMethod))(iSample), results.(append(iMethod, 'Time'))(iSample)] = getPCDR3(s1, s2, mx, Sigmax, iMethod);
    end
    
    s1Array(:,iSample) = [s1.a, s1.q, s1.tc', s1.eps];
    s2Array(:,iSample) = [s2.a, s2.q, s2.tc', s2.eps];
end
close(f)

%sort data by PCD-exact in ascent way and store sorted data in a struct
%array format
[results.Exact, index] = sort(results.Exact(1,1:sampleNumber), 2);

for i = 2:size(keySet,2)
    iMethod = string(keySet(i));
    results.(genvarname(iMethod)) = results.(genvarname(iMethod))(index);
    results.(genvarname(append(iMethod, 'Time'))) = results.(genvarname(append(iMethod, 'Time')))(index);
end

% s1 s2 Information
results.s1Array = s1Array(:, index);
results.s2Array = s2Array(:, index);
end