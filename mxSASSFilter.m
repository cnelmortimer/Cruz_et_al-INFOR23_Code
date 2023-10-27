function [filteredDB, evalsConsumed] = mxSASSFilter(DBPath, database, matMolA, weightsA, ovA, nDims, params)
    qnt = params(1);
    localEvals = params(2);
    
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
            func = @(x) - mxObjFunc(x, matMolA, weightsA, ovA, matMolB, weightsB, ovB, boundsAB); % We are maximizing in reality
            [~, localRank] = MySASS(filterX, - bufferVal(j), func, localEvals-1);
            bufferVal(j) = -localRank;
        end
        preVal(i) = max(bufferVal); % Considering the best insight
    end
    
    if qnt > 1 % qnt is a direct quantity of the most promising (best averaged) molecules
        [~, indices] = maxk(preVal, qnt);
    else % qnt is between 0 and 1, and refers to a degradation percentage related to the best preliminary result
        indices = find(preVal >= ( (1-qnt) * max(preVal) ) );
    end
    
    filteredDB = database(indices);
    evalsConsumed = numFilters*localEvals*numel(database);
end

% ------ Internal auxiliary functions: ------

function filters = loadFilters(nDims)
    if nDims==10
        filters = [0 1 0 0 0 0 0 0 0 0; ...     % noMod
                pi 1 0 0 0 0 0 0 0 0; ...    % x180
                pi 0 1 0 0 0 0 0 0 0; ...    % y180
                pi 0 0 1 0 0 0 0 0 0];       % z180
    elseif nDims==6
	filters = [0 0 pi/2 0 0 0; ...  % No mod
                pi 0 pi/2 0 0 0; ...    % Around X axis
                pi pi/2 pi/2 0 0 0; ... % Around Y axis
                pi 0 0 0 0 0];          % Around Z axis
    else
       error('The filters for the 6-vars version are knot known yet'); 
    end
end

function [currentPoint, currentVal] = MySASS(normal_startPoint, startValue, normal_func, MaxEvals)
    siz = size(normal_startPoint);

    MaxSigma = 1.0;
    MinSigma = 1e-5;
    sigma = MaxSigma;
    Scnt = 5; % Consecutive success trigger
    Fcnt = 3; % Consecutive fail trigger
    ex = 2.0; % Expansion factor
    ct = 0.5; % Contraction factor
    
    currentPoint = normal_startPoint; % Move if and only if it improves
    currentVal = startValue;
    
    bias = zeros(siz); 
    
    scnt = 0; % Number of consecutive successes
    fcnt = 0; % Number of consecutive failures
    evals = 0; % Number of times that 'normal_func' has been called
    
    while(evals<MaxEvals)
        perturbation = bias + normrnd(0, sigma, siz); % Normal Perturbation (N(0, sigma) + bias == N(bias, sigma)). <With size of bias>
                                                          % Ex: rng(4321); normrnd(0, 1) + 7; ans = 6.2424; rng(4321); normrnd(7, 1); ans = 6.2424
                                                          
        x_prime = currentPoint + perturbation; % POSITIVE DIRECTION
        
        x_prime(x_prime < 0) = 0; % Fixing the lower bounds
        x_prime(x_prime > 1) = 1; % Fixing the upper bounds
        
        val_prime = normal_func(x_prime);
        evals = evals + 1; % One more function evaluation
        
        if val_prime < currentVal
            currentPoint = x_prime; % Move to the improving point!
            currentVal = val_prime;
            bias = 0.2*bias + 0.4*perturbation; % perturbation = dif in N(0, sigma) + bias = dif in N(bias, sigma)
            scnt = scnt + 1; % One more success
            fcnt = 0;        % Reset the number of failures
        elseif(evals < MaxEvals)
            x_prime = currentPoint - perturbation; % NEGATIVE DIRECTION
            
            x_prime(x_prime < 0) = 0; % Fixing the lower bounds
            x_prime(x_prime > 1) = 1; % Fixing the upper bounds
        
            val_prime = normal_func(x_prime);
            evals = evals + 1; % One more function evaluation
            
            if val_prime < currentVal
                currentPoint = x_prime; % Move to the improving point!
                currentVal = val_prime;
                bias = bias - 0.4*perturbation;
                scnt = scnt + 1; % One more success
                fcnt = 0;        % Reset the number of failures
            else
                bias = 0.5*bias;
                fcnt = fcnt + 1; % One more fail
                scnt = 0;        % Reset the number of successes
            end
        end
        
        if scnt > Scnt
            sigma = ex*sigma; % Expansion % scnt = 0; % Do not reset the number of successes to allow premature exit
        elseif fcnt > Fcnt
            sigma = ct*sigma; % Contraction % fcnt = 0; % Do not reset the number of failures to allow premature exit
        end
        if sigma < MinSigma || fcnt>20 % Restart the search
            sigma = MaxSigma;
            scnt = 0; 
            fcnt = 0; 
            bias = 0*bias;
        elseif sigma > MaxSigma % Do not move too fast
            sigma = MaxSigma;        
        end
    end
end
