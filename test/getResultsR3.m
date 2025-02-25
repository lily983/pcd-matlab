function results = getResultsR3(sampleNumber, objectType, Sigma1, Sigma2, keySet)
% getResultsR3 Generate random objects s1 and s2 and return PCD value
% generate by different methods
%
% Inputs
%     sampleNumber      number of random samples
%     objectType              object is sphere or ellipsoid or superquadrics
%     Sigma1  position error covariance matrix of object 1
%     Sigma2  position error covariance matrix of object 2
%     keySet   list of methods names
%
% Outputs
%     results                           struct array variables, including
%     collision status and PCD results by each methos

f = waitbar(0, 'Starting');

% Check if two objects 
if all(diag(Sigma1))
    TwoErrorCase = 1;
    filePath = [objectType  '-TwoError-'];
else
    TwoErrorCase = 0;
    filePath = [objectType  '-SingleError-'];
end

%  Create .jason file
currentDateTime = now;
formattedDateTime = datestr(currentDateTime, 'yyyy-mm-dd HH:MM:SS');
filePath = [filePath  formattedDateTime  '.json'];

iSample = 1;
while iSample <= sampleNumber
    waitbar(iSample/sampleNumber, f, sprintf('Progress: %d %%', floor(iSample/sampleNumber*100)));
    
    if strcmp(objectType, 'sphere')
        SN = [10, 10];
        
        s1 = SuperQuadrics({(rand(1)+0.2)*ones(1,3), [1,1], [0, 0]...
            rand(3,1).*0.1, getRandomQuaternion(), SN});
        
        s2 = SuperQuadrics({(rand(1)+0.2)*ones(1,3),  [1,1], [0, 0]...
            rand(3,1)+0.3, getRandomQuaternion(),SN});
        
    elseif strcmp(objectType, 'ellipsoid')
        SN = [10, 10];
        
        s1 = SuperQuadrics({rand(1,3)+0.2, [1,1], [0, 0]...
            rand(3,1).*0.1, getRandomQuaternion(), SN});
        
        s2 = SuperQuadrics({rand(1,3)+0.2,  [1,1], [0, 0]...
            rand(3,1)+0.3, getRandomQuaternion(), SN});
        
    elseif strcmp(objectType, 'superquadrics')
        SN = [40, 40];
        
        s1 = SuperQuadrics({rand(1,3)+0.2, [(rand(1)+0.01)*2, rand(1)*2], [0, 0]...
            rand(3,1).*0.1, getRandomQuaternion(), SN});
        
        s2 = SuperQuadrics({rand(1,3)+0.2, [(rand(1)+0.01)*2, rand(1)*2], [0, 0]...
            rand(3,1)+0.3, getRandomQuaternion(), SN});
    end
    
    try 
        [flag, dist, ~, condition] = collision_cfc(s1,s2, 'least-squares');
    catch 
        iSample = iSample-1;
        continue;
    end
    
    % Check cfc condition between s1 and s2, if not converge, then generate
    % new pairs
    if objectType=="superquadrics"
        if condition > 1e-02 || isnan(condition)
            if iSample ==1
                iSample = 1;
            end
            iSample = iSample-1;
            continue;
        end
    end
    
    flagArray(iSample) = round(flag);
    distArray(iSample) = dist;
    
    % The mean of the relative position error
    mx = s2.tc-s1.tc;
    
    % Now transform the position error covariance matrix from the
    % body-fixed frame to the world frame
    R1 = quat2rotm(s1.q);
    R2 = quat2rotm(s2.q);
    % ellipsoid: x' * inv(Sigma) * x =1 => x' * inv(R*diag^2*R') * x=1
    Sigma1 = R1 * Sigma1 * R1';
    Sigma2 = R2 * Sigma2 * R2';
    
%     % Check if two objects 
%     if all(diag(Sigma1))
%         TwoErrorCase = 1;
%     else
%         TwoErrorCase = 0;
%     end
    
    for i = 1:size(keySet,2)
        iMethod = string(keySet(i));
        
%         if two object collide at mean pose, the following algorithms will
%         not work
        if flag && any(strcmp(iMethod, {'Maxpdf', 'Maxpdf_SQ', 'LCC_center_point', 'LCC_closed_point', ...
                'LCC_tangent_point', 'Quadratic_exact', 'Quadratic_bound', 'Divergence_mesh', 'LCC_center_point_cfc', 'LCC_tangent_point_cfc'}))
            results.(genvarname(iMethod))(iSample)=1;
            results.(append(iMethod, 'Time'))(iSample)=NaN;
            continue;
        end
        
        % If TwoErrorCase, only get results for following methods
        if TwoErrorCase 
            if any(strcmp(iMethod, {'Fast_sampling', 'Exact_two', 'EB_99', 'EB_95', 'GMF', 'Quadratic_exact', 'Quadratic_bound',...
                    'LCC_center_point', 'LCC_closed_point', 'LCC_tangent_point', 'LCC_center_point_cfc', 'LCC_tangent_point_cfc', 'Divergence_mesh'})) 
                [results.(genvarname(iMethod))(iSample), results.(append(iMethod, 'Time'))(iSample)] ...
                    = getPCDR3(s1, s2, mx, Sigma1, Sigma2, iMethod);
            else
                results.(genvarname(iMethod))(iSample)=NaN;
                results.(append(iMethod, 'Time'))(iSample)=NaN;
            end
        % If not TwoErrorCase, the following methods will not work, omit
        % them
        else
            if any(strcmp(iMethod, {'Fast_sampling','Exact_two'})) 
                results.(genvarname(iMethod))(iSample)=NaN;
                results.(append(iMethod, 'Time'))(iSample)=NaN;
            else
                [results.(genvarname(iMethod))(iSample), results.(append(iMethod, 'Time'))(iSample)] ...
                    = getPCDR3(s1, s2, mx, Sigma1, Sigma2, iMethod);
            end
        end
    end
    
    s1Array(:,iSample) = [s1.a, s1.q, s1.tc', s1.eps];
    s2Array(:,iSample) = [s2.a, s2.q, s2.tc', s2.eps];
    
    % Write results
    jsonData = jsonencode(results);
    fileID = fopen(filePath, 'w');
    fprintf(fileID, '%s', jsonData);
    
    iSample = iSample + 1;
end
close(f)

% Close the file
fclose(fileID);

%sort data by Exact/Exact_two in ascent way and store sorted data in a struct
%array format
if TwoErrorCase
%     baselineName = 'Exact_two';
    baselineName = 'Fast_sampling';
    [results.(genvarname(baselineName)), index] = sort(results.(genvarname(baselineName)), 2);
else
    baselineName = 'Exact';
    [results.(genvarname(baselineName)), index] = sort(results.(genvarname(baselineName)), 2);
end

for i = 1:size(keySet,2)
    iMethod = string(keySet(i));
    if strcmp(iMethod, baselineName)
        continue
    end
    results.(genvarname(iMethod)) = results.(genvarname(iMethod))(index);
    results.(genvarname(append(iMethod, 'Time'))) = results.(genvarname(append(iMethod, 'Time')))(index);
end

% s1 s2 Information
results.s1Array = s1Array(:, index);
results.s2Array = s2Array(:, index);

% flag and dist
results.flag = flagArray(:, index);
results.dist = distArray(:, index);
end