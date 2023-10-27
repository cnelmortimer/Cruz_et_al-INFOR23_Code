function output = mxConfParDBExplorer(DBPath, refMolName, nDims, optHandle, optParams, histoLength, filterHandle, filterParams)
    database = dir([DBPath, '/DB*.mol2']);
    if isempty(database)
       error('The database cannot be read');
    end

    [database, refNames] = getRefs(database, refMolName);
    if isempty(refNames)
        error('No reference files were found');
    end

    pool = gcp();
    numWorkers = pool.NumWorkers;

    for r=1:1:numel(refNames)
        currFileRef = refNames{r};
        currFileName = erase(currFileRef, '.mol2');
        % Let us fix the reference:
        matMol1 = mxReadMol([DBPath, '/', currFileRef]);
        weights1 = mxCalculateWeights(matMol1);
        ov1 = mxPreciseOverlapWEGA(matMol1, weights1, matMol1, weights1);

        woHistory = cell(numWorkers, 3);% A place per worker
        
        partialDB = database;
        if exist('filterHandle', 'var') && exist('filterParams', 'var')
            partialDB = filterHandle(DBPath, database, matMol1, weights1, ov1, nDims, filterParams);
        end

        parfor id=1:numWorkers
            woHistory(id, :) = { repmat({'0'}, 1, histoLength), repmat({zeros(1, 10)}, 1, histoLength), -inf*ones(1, histoLength) };
            for j=id:numWorkers:numel(partialDB)
                matMol2 = mxReadMol([DBPath, '/', partialDB(j).name]);
                weights2 = mxCalculateWeights(matMol2);
                ov2 = mxPreciseOverlapWEGA(matMol2, weights2, matMol2, weights2);
                bounds = mxComputeBounds(matMol1, matMol2, nDims);
                func = @(x) - mxObjFunc(x, matMol1, weights1, ov1, matMol2, weights2, ov2, bounds); % We are maximizing in reality

                [x, fval] = optHandle(nDims, func, optParams{:});
                woHistory(id, :) = compareAndSave(woHistory(id, :), histoLength, erase(partialDB(j).name, '.mol2'), x, -fval);
            end
        end
        [foundMol, foundX, foundVal] = mixHistories(woHistory, histoLength);
        output.(currFileName) = {foundMol, foundX, foundVal};
    end
end

% ----------------------- Internal auxiliary functions

function [database, refNames] = getRefs(database, refMolName)
    refNames = cell(0);
    moved = [];
    for i=1:1:numel(database)
        if contains(database(i).name, refMolName)
            refNames(end + 1) = {database(i).name};
            moved = [moved, i];
        end
    end
    database(moved) = [];
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
