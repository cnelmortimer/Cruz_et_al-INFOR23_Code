function [filteredDB, evalsConsumed] = mxDBFilter(DBPath, database, matMolA, weightsA, ovA, nDims, qnt)
    filters = loadFilters(nDims);
    numFilters = size(filters, 1); % As many filters as rows in 'filters'
    
    preVal = zeros(numel(database), 1);
    
    parfor i=1:numel(database) % Remember not to keep the reference in the DB!
        matMolB = mxReadMol([DBPath, '/', database(i).name]);
        weightsB = mxCalculateWeights(matMolB);
        ovB = mxPreciseOverlapWEGA(matMolB, weightsB, matMolB, weightsB);
        
        boundsAB = mxComputeBounds(matMolA, matMolB, nDims);
        denom = boundsAB(2,:) - boundsAB(1,:); % max - min
        probDims = (denom==0); % Problematic dimensions => Lower bound = Upper bounds = No diff
        bufferVal = zeros(1, numFilters);
        for j=1:1:numFilters
            filterX = (filters(j,:) - boundsAB(1,:)) ./ denom; % Normalize this filter for this context
            if any(probDims)
                filterX(probDims) = boundsAB(1, probDims); % Avoid no-diff NaN's by setting them to their lower bound
            end
            bufferVal(j) = mxObjFunc(filterX, matMolA, weightsA, ovA, matMolB, weightsB, ovB, boundsAB);
        end
        preVal(i) = max(bufferVal); % Considering the best insight
    end
    
    if qnt > 1 % qnt is a direct quantity of the most promising (best averaged) molecules
        [~, indices] = maxk(preVal, qnt);
    else % qnt is between 0 and 1, and refers to a degradation percentage related to the best preliminary result
        indices = find(preVal >= ( (1-qnt) * max(preVal) ) );
    end
    
    filteredDB = database(indices);
    evalsConsumed = numFilters*numel(database);
end

% ------ Internal auxiliary functions: ------

function filters = loadFilters(nDims)
    if(nDims==10)
        filters = [0 1 0 0 0 0 0 0 0 0; ...   % noMod
                   pi 1 0 0 0 0 0 0 0 0; ...  % x180
                   pi 0 1 0 0 0 0 0 0 0; ...  % y180
                   pi 0 0 1 0 0 0 0 0 0];     % z180
    elseif(nDims==6)
        filters = [0 0 pi/2 0 0 0; ...      	% No mod
                pi 0 pi/2 0 0 0; ...     	% Around X axis
                pi pi/2 pi/2 0 0 0; ... 	% Around Y axis
                pi 0 0 0 0 0];            	% Around Z axis
    else
       error('Unexpected number of variables!'); 
    end
end
