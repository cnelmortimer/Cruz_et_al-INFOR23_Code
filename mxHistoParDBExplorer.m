function [foundMol, foundX, foundVal] = mxHistoParDBExplorer(DBPath, refMolName, nDims, optHandle, optParams, histoLength, filterHandle, filterParams)
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
    
    pool = gcp(); % parpool();
    numWorkers = pool.NumWorkers;
    woHistory = cell(numWorkers, 3);% A place per worker
    
    parfor id=1:numWorkers
        woHistory(id, :) = { repmat({'0'}, 1, histoLength), repmat({zeros(1, 10)}, 1, histoLength), -inf*ones(1, histoLength) };
        for j=id:numWorkers:numel(database)
            matMol2 = mxReadMol([DBPath, '/', database(j).name]);
            weights2 = mxCalculateWeights(matMol2);
            ov2 = mxPreciseOverlapWEGA(matMol2, weights2, matMol2, weights2);
            bounds = mxComputeBounds(matMol1, matMol2, nDims);            
            func = @(x) - mxObjFunc(x, matMol1, weights1, ov1, matMol2, weights2, ov2, bounds); % We are maximizing in reality
                
            [x, fval] = optHandle(nDims, func, optParams{:});
            woHistory(id, :) = compareAndSave(woHistory(id, :), histoLength, erase(database(j).name, '.mol2'), x, -fval);
        end
    end
        
    [foundMol, foundX, foundVal] = mixHistories(woHistory, histoLength);
    % delete(pool);
end

% ----------------------- Internal auxiliary functions

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

function hist = compareAndSave(hist, histLength, newName, newX, newVal)
    names = hist{1};
    xs = hist{2};
    vals = hist{3};
    
    focus = find(newVal>=vals, 1); % Preliminary shortcut
    if isempty(focus)
       return; % Go out without modifying the history
    else % Inject:
        for i=histLength:-1:focus+1
            vals(i) = vals(i-1);
            names(i) = names(i-1);
            xs(i) = xs(i-1);
        end
        vals(focus) = newVal;
        names(focus) = {newName};
        xs(focus) = {newX};
    end
    hist = {names, xs, vals};
end

function [bestNames, bestXs, bestVals] = mixHistories(history, histoLength)    
    bestNames = horzcat(history{:, 1})';
    bestXs = horzcat(history{:, 2})';
    bestVals = horzcat(history{:, 3})';
    
    [~, indices] = sort(bestVals, 'descend');
    indices = indices(1:histoLength); % As every thread has the same histoLength, we must limit to one time this quantity
    bestVals = bestVals(indices);
    bestNames = bestNames(indices);
    bestXs = cell2mat(bestXs(indices));
end
