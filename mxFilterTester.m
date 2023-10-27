function mxFilterTester(filterHandle, nDims, filterParams)
    [DBPath, fullDB, targets] = loadContext();

    safetyChecks = zeros(numel(targets), 1);
    accumMols = 0;
    accumEvals = 0;
    for i=1:1:numel(targets)
        refMolName = targets{i}{1};
        desireName = targets{i}{2};
        disp(['Working with: ', refMolName]);
        
        matMol1 = mxReadMol([DBPath, '/', refMolName, '.mol2']);
        weights1 = mxCalculateWeights(matMol1);
        ov1 = mxPreciseOverlapWEGA(matMol1, weights1, matMol1, weights1);
        
        partialDB = selfExclude(fullDB, refMolName);
        numTotal = numel(partialDB);
        
        [filteredDB, evalsConsumed] = filterHandle(DBPath, partialDB, matMol1, weights1, ov1, nDims, filterParams);
        numKept = numel(filteredDB);
        accumMols = accumMols + numKept;
        accumEvals = accumEvals + evalsConsumed; % (As the numer may differ when dealing with conformations)
        
        if checkSuccess(filteredDB, desireName)
            safetyChecks(i) = 1;
        else
            disp('Warning: Goal lost :-/');
        end
        disp(['Kept: ', num2str(numKept), ' (out of ', num2str(numTotal), ') => Skipped: ', num2str(numTotal - numKept), '. Compression rate: ', num2str( (numTotal - numKept) / numTotal)]);
        disp('------------------------');
    end
    disp(['END. Total evaluations: ', num2str(accumEvals), '. Average number kept: ', num2str(accumMols/numel(targets)), '. Succes: ', num2str(sum(safetyChecks)), ' (', num2str( sum(safetyChecks) / numel(targets) ), ')']);
end

% ----------------------- Internal auxiliary functions

function [DBPath, fullDB, targets] = loadContext()
    DBPath = '../INPUT/miriamfda';
    targets = { {'DB00529', 'DB09294'}, {'DB00331', 'DB09210'}, {'DB01352', 'DB00306'}, {'DB01365', 'DB00191'}, {'DB00380', 'DB01041'}, ...
    {'DB06216', 'DB00370'}, {'DB00693', 'DB01619'}, {'DB07615', 'DB00721'}, {'DB09219', 'DB01320'}, {'DB00674', 'DB01619'}, ...
    {'DB01198', 'DB00402'}, {'DB00887', 'DB00837'}, {'DB00246', 'DB01261'}, {'DB00381', 'DB01023'}, {'DB09237', 'DB01054'}, ...
    {'DB00876', 'DB09039'}, {'DB00254', 'DB00595'}, {'DB00351', 'DB04839'}, {'DB01196', 'DB00286'}, {'DB01621', 'DB01148'}, ...
    {'DB09236', 'DB01054'}, {'DB08903', 'DB00333'}, {'DB00632', 'DB00464'}, {'DB01419', 'DB06605'}, {'DB00320', 'DB00728'}, ...
    {'DB00728', 'DB01339'}, {'DB00503', 'DB00701'}, {'DB01232', 'DB00212'}, {'DB00309', 'DB00541'}, {'DB04786', 'DB00511'}, ...
    {'DB09114', 'DB08993'}, {'DB06439', 'DB00207'}, {'DB01078', 'DB00511'}, {'DB01590', 'DB00877'}, {'DB04894', 'DB00646'}, ...
    {'DB00403', 'DB08874'}, {'DB00732', 'DB06287'}, {'DB00050', 'DB00569'}, {'DB06699', 'DB09099'}, {'DB06219', 'DB00512'} };

    fullDB = dir([DBPath, '/DB*.mol2']);
    if isempty(fullDB)
       error('The database cannot be read'); 
    end
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

function success = checkSuccess(namesDB, desiredName)
    keptNames = vertcat(namesDB.name);
    numNames = size(keptNames, 1); % As many names as rows
    i = 1;
    success = false;
    while i<=numNames
        if contains(keptNames(i, :), desiredName) % Using contains in case we work with conformations, i.e., base name plus variations
            success = true;
            break;
        end
        i = i + 1;
    end
end
