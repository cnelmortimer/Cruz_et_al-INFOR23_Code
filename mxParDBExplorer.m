function [foundMol, foundX, foundVal] = mxParDBExplorer(DBPath, refMolName, nDims, optHandle, optParams, filterHandle, filterParams)
    database = dir([DBPath, '/DB*.mol2']);
    if isempty(database)
       error('The database cannot be read'); 
    end
        
    % Let us fix the reference:
    matMol1 = mxReadMol([DBPath, '/', refMolName, '.mol2']);
    weights1 = mxCalculateWeights(matMol1);
    ov1 = mxPreciseOverlapWEGA(matMol1, weights1, matMol1, weights1);
    
    database = selfExclude(database, refMolName); % Do not compare to itself!
    if exist('filterHandle', 'var') && exist('filterParams', 'var')
        database = filterHandle(DBPath, database, matMol1, weights1, ov1, nDims, filterParams);
    end
    
    pool = parpool();
    numWorkers = pool.NumWorkers;
    woFoundMol = repmat({'0'}, 1, numWorkers);% A place per worker
    woFoundX = repmat({zeros(1, 10)}, 1, numWorkers);
    woFoundVal = -inf*ones(1, numWorkers);
    
    parfor id=1:numWorkers
        for j=id:numWorkers:numel(database)            
            matMol2 = mxReadMol([DBPath, '/', database(j).name]);
            weights2 = mxCalculateWeights(matMol2);
            ov2 = mxPreciseOverlapWEGA(matMol2, weights2, matMol2, weights2);
            bounds = mxComputeBounds(matMol1, matMol2, nDims);            
            func = @(x) - mxObjFunc(x, matMol1, weights1, ov1, matMol2, weights2, ov2, bounds); % We are maximizing in reality
                
            [x, fval] = optHandle(nDims, func, optParams{:});
            if -fval > woFoundVal(id)
                woFoundVal(id) = -fval;
                woFoundX{id} = x;
                woFoundMol{id} = erase(database(j).name, '.mol2');
            end
        end
    end
        
    [foundVal, idx] = max(woFoundVal);
    foundMol = woFoundMol{idx};
    foundX = woFoundX{idx};
    delete(pool);
end

function database = selfExclude(database, refMolName)
    i = 1;
    while i <= numel(database)
        if contains(database(i).name, refMolName)
            database(i) = []; % Remove this item;
        else
            i = i + 1;
        end
    end
end
