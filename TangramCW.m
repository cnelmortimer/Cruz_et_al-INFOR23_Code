function [x,fval] = TangramCW(nDim, normalFunc, maxEvals, localBudget, wrappedVars)
%TANGRAM A derivative-free stochastic optimizer with few parameters
% <EXPERIMENTAL VERSION using centroids instead of midpoints>
% Tangram is a derivative-free stochastic optimizer with few parameters:
% nDim => Number of dimensions of the search space (and params for normalFunc)
% normalFunc => Objective function assuming it expects all the params in [0, 1] (normalized)
% maxEvals => Maximum number of function evaluations
% localBudget => Number of function evaluations allowed to every SASS run
% wrappedVars => Indices of the variables to wrap around instead of truncating
% It returns the best (normalized) point found and its value according to normalFunc
    availableEvals = maxEvals;
    if ~exist('localBudget', 'var') || isempty(localBudget)
        localBudget = 32;
    end
    
    [facets, numFacets] = getFacets(nDim); % Getting the facets used for computing every centroid
    domainDiameter = sqrt(nDim); % sqrt( (1_1 - 0_1)^2 + (1_2 - 0_2)^2 + ... (1_nDim - 0_nDim)^2) )
    
    % Step 1: Creating the initial point:
    x = 0.5*ones(1, nDim);
    fval = normalFunc(x);
    availableEvals = availableEvals - 1;
    
    incisive = false;
    if maxEvals < (localBudget + 1 + numFacets) % In general: maxEvals < 33 + 2*N
        incisive = true; % Incisive        
    end
    
    newPoints = zeros(numFacets, nDim); % Memory room reserved once!
    newVals = zeros(numFacets, 1) + inf;
    newRadius = zeros(numFacets, 1);
    
    while(availableEvals > 0)
        % Step 2: Global phase -> allowing SASS to take long steps
        budget = min(availableEvals, localBudget);
        
        [x, fval] = MySASS(x, fval, normalFunc, budget, domainDiameter, wrappedVars); % Move if and only if it improves (bound-checking integrated)
        availableEvals = availableEvals - budget;
        
        % Step 3: (mode 1 (Standard) => Division) || (mode 2 (Incisive) => Division & Local Phase)
        newPoints = 0*newPoints; % Resetting for this iteration
        newVals = newVals + inf;
        newRadius = 0*newRadius;
        if(~incisive)
            [newPoints, newVals, newRadius, availableEvals] = StandardDivisionStage(x, facets, numFacets, newPoints, newVals, newRadius, availableEvals, normalFunc);
            [newPoints, newVals, availableEvals] = StandardLocalSearch(newPoints, newVals, newRadius, availableEvals, localBudget, normalFunc, wrappedVars); % Step 4: Local search
        else
            [newPoints, newVals, availableEvals] = IncisiveDivisionAndSearch(x, facets, numFacets, newPoints, newVals, newRadius, availableEvals, localBudget, normalFunc, wrappedVars);
        end
        
        % Step 5 (4 in Incisive): Update the center:
        [candidateMin, minPosition] = min(newVals);
        if candidateMin < fval
            fval = candidateMin;
            x = newPoints(minPosition, :);
        end
    end
end

%-------------------Internal auxiliary functions:--------------------------

function [facets, numFacets] = getFacets(nDim)
    corners = ff2n(nDim); % Fast way to get all the binary combinations from 0 to 2^N - 1
    numFacets = 2*nDim;
    for i=1:1:nDim
       indices = (corners(:, i)==1);
       facets{2*i - 1} = corners(indices, :); %#ok<AGROW>
       indices = (corners(:, i)==0);
       facets{2*i} = corners(indices, :); %#ok<AGROW>
    end
end

function [newPoints, newVals, newRadius, availableEvals] = StandardDivisionStage(x, facets, numFacets, newPoints, newVals, newRadius, availableEvals, normalFunc)
    k = size(facets{1}, 1) + 1; % Every centroid refers to k points = The corners involved & x = (2^(nDim-1)) + 1
    for i=1:1:numFacets
       newPoints(i, :) = (sum(facets{i}) + x) / k; % Getting the centroid of the reference & the facets
       newRadius(i) = norm(newPoints(i, :) - x);
       if availableEvals > 0
            newVals(i) = normalFunc(newPoints(i,:));
            availableEvals = availableEvals - 1;
        else
            break;
        end
    end
end

function [newPoints, newVals, availableEvals] = StandardLocalSearch(newPoints, newVals, newRadius, availableEvals, localBudget, normalFunc, wrappedVars)
    [~, indices] = sort(newVals, 'ascend'); % Better (lower) regions are explored first
    for i=indices'
        budget = min(availableEvals, localBudget);
        if budget<=0
            break;
        else
            [newPoints(i, :), newVals(i)] = MySASS(newPoints(i, :), newVals(i), normalFunc, budget, newRadius(i), wrappedVars);
            availableEvals = availableEvals - budget;
        end
    end
end

function [newPoints, newVals, availableEvals] = IncisiveDivisionAndSearch(x, facets, numFacets, newPoints, newVals, newRadius, availableEvals, localBudget, normalFunc, wrappedVars)
    k = size(facets{1}, 1) + 1; % Every centroid refers to k points = The corners involved & x = (2^(nDim-1)) + 1
    for i = 1:1:numFacets
        newPoints(i, :) = (sum(facets{i}) + x) / k;
        newRadius(i) = norm(newPoints(i, :) - x);
        if availableEvals > 0
            newVals(i) = normalFunc(newPoints(i,:));
            availableEvals = availableEvals - 1;
                
            budget = min(availableEvals, localBudget);
            if budget<=0
                break;
            else
                [newPoints(i, :), newVals(i)] = MySASS(newPoints(i, :), newVals(i), normalFunc, budget, newRadius(i), wrappedVars);
                availableEvals = availableEvals - budget;
            end
        else
            break;
        end
    end
end

function [currentPoint, currentVal] = MySASS(normal_startPoint, startValue, normal_func, MaxEvals, MaxStepSize, wrappedVars)
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
        amplitude = norm(perturbation);
        if amplitude > MaxStepSize % norm( (currentPoint + perturbation) - currentPoint) = norm(perturbation)
            perturbation = MaxStepSize*(perturbation/amplitude);
        end
        
        x_prime = currentPoint + perturbation; % POSITIVE DIRECTION
        x_prime = fixLimits(x_prime, wrappedVars);
        
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
            x_prime = fixLimits(x_prime, wrappedVars);
        
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

function vec = fixLimits(vec, wrappedVars)
    for i = wrappedVars
        if vec(i)>1 || vec(i)<0
            part = abs(vec(i)) - floor(abs(vec(i))); % Taking the decimal part (https://www.mathworks.com/matlabcentral/answers/21339-how-to-seperate-fractional-and-decimal-part-in-a-real-number)
            if vec(i)< 0
                vec(i) = 1 - part; % Starting from 1
            else
                vec(i) = part; % Starting from 0
            end
        end
    end
    vec(vec < 0) = 0; % The previous (wrapped) parts should not be affected now
    vec(vec > 1) = 1; % Fixing the bounds
end
