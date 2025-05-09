function [ table_outputs, table_ee, factors_over, n, points, initialization ] = moat_loop_function( nfac, studied_function, table_outputs, table_ee, factors_over, n, points, initialization, varargin )
% Function for the calculation of the algorithm at one step.
% This function must be used in a loop (as in the file program_to_run.m).
% See program_to_run.m for further information about the variables.
%
% 8 INPUTS
% 1) nfac : Number of factors of uncertainty.         
% 2) studied_function : Function studied.                 
% 3) table_outputs : All the outputs of the simulations runs.
% 4) table_ee : All the elementary effects.
% 5) factors_over : Indexes of the factors over the limit at the end of the last step.
% 6) n : Step number.
% 7) points : Sampled points of the factors hyperspace.         
% 8) initialization : Boolean, the calculations foccus on the factors under the limit if it equals 1.
%
% 6 OUTPUTS : table_outputs2, table_ee2, factors_over2, n2, points2, initialization2 (see the description of the corresponding inputs)

if nargin>8
    options = varargin{1};
else
    options.limit_factor = 0.1; % default limit factor
    options.initialization_factor = 0.5; % default factor for initialization
end

%% DENSIFICATION
% Information about the points sampled
npoint = size(points,1); % npoint is the number of points already sampled.
% If the step number is larger than the number of points already sampled, new points must be sampled
if npoint < n
    if npoint==0 % At the beginning of the algorithm, 4 points are first sampled
        npoint = 4;
    else
        npoint = npoint*3; % The number of points is multiplied by 3 at each densification
    end
    
    delta=1/npoint; % Distance between two successive discrete values in a given dimension.
    mini=delta/2; % Minimum discrete value in a given dimension.
    
    % Values of the npoint discrete values in each dimension
    coord = (0.5*delta):delta:1;

    if npoint>4 % If points have been sampled before
        % Values which have already been used
        npoint_prec=npoint/3;
        delta_prec=1/npoint_prec;
        mini_prec=delta_prec/2;
        maxi_prec=mini_prec+delta_prec*(npoint_prec-1);
        coord_prec = 0:npoint_prec-1;
        coord_prec = coord_prec/(npoint_prec-1)*(maxi_prec-mini_prec)+mini_prec;
        
        puissance_verifiee = -round(log10(mini)-6); % Accuracy used to compare values.
        
        [~,ia] = setdiff(round(coord*10^puissance_verifiee),round(coord_prec*10^puissance_verifiee)); % do not repeat previous calculations
        coord = coord(ia);
    end
    
    % Sampling the new points in a Latin hypercube : For each factor (for a given column), each of the new coordinates is one different value of coord
    for i=nfac:-1:1
        points_extra(:,i) = coord(randperm(length(coord))); % Additional points
    end
    
    points = [points; points_extra]; % The additional points are added to the former points
end
%% CALCULATION OF THE ELEMENTARY EFFECTS
% Variation applied to the factors to calculate the elementary effects
if mod(n,2) == 1 % If the step number is odd
    variation = 0.5;
else % If the step number is even
    variation = 0.75*0.5; 
end
lines_outputs = size(table_outputs,1); % lines_outputs is the maximum step number reached
if ~initialization % At the beginning, the elementary effects are calculated for all the factors
    table_outputs(n,1) = studied_function(points(n,:)); % Output of the function at the n-th sampled point.
    
    for i=1:nfac % For each factor
        if points(n,i) < 0.5 % If its coordinate if smaller than 0.5, a positive variation is applied
            table_outputs(n,1+i) = studied_function([points(n,1:i-1) points(n,i)+variation points(n,i+1:nfac)]);
        else % If its coordinate if larger than 0.5, a negative variation is applied
            table_outputs(n,1+i) = studied_function([points(n,1:i-1) points(n,i)-variation points(n,i+1:nfac)]);
        end
        table_ee(n,i) = abs(table_outputs(n,1+i)-table_outputs(n,1))/variation; % Elementary effect of the i-th factor at the n-th sampled point.
    end
    
    if n>=3 % From the third step, the ability to define a limit between the factors with a low and a high influence is estimated
        test_initialization = zeros(3,nfac);
        for i=1:3 % For three steps
            test_initialization(i,:) = max(table_ee(1:n-3+i,:),[],1); % Maxima of the elementary effects obtained for each factor (column) at each of the last three steps (three lines)
        end
        test_initialization_sort = sort(test_initialization,2); % Ordering the maxima of the elementary effects of each factor.
        test_initialization_delta_max_sort = sort(max(diff(test_initialization_sort,1,2),[],2),'descend'); % sorted max Variation between the successive maxima of each step
        if (test_initialization_delta_max_sort(3)>=options.initialization_factor*test_initialization_delta_max_sort(1)) % If the smallest value is larger than half of the largest value, the variations are considered to be roughly stable.
           initialization = true; % The calculations will then foccus on the factors under the limit defined with the largest variation between successice maxima.
        end
    end
    
elseif n<=lines_outputs % If the algorithm returned to a previous step
    
    for i=1:nfac % For each of the factors
        if (table_ee(n,i)==-1) && ~ismember(i,factors_over) % If the elementary effect has not been calculated before and if the factor is not over the limit         
            if points(n,i) < 0.5 % If its coordinate if smaller than 0.5, a positive variation is applied
                table_outputs(n,1+i) = studied_function([points(n,1:i-1) points(n,i)+variation points(n,i+1:nfac)]);
            else % If its coordinate if larger than 0.5, a negative variation is applied
                table_outputs(n,1+i) = studied_function([points(n,1:i-1) points(n,i)-variation points(n,i+1:nfac)]);
            end
            table_ee(n,i) = abs(table_outputs(n,1+i)-table_outputs(n,1))/variation; % Elementary effect of the i-th factor at the n-th sampled point.
        end
    end
    
else % If it is a new step where the calculations have to be foccused on the factors under the limit
    table_outputs(n,1) = studied_function(points(n,:)); % Output of the function at the n-th sampled point.
    for i=1:nfac % For each factor
        if ~ismember(i,factors_over) % If the factor is not over the limit   
            if points(n,i) < 0.5 % If the elementary effect has not been calculated before and if the factor is not over the limit
                table_outputs(n,1+i) = studied_function([points(n,1:i-1) points(n,i)+variation points(n,i+1:nfac)]);
            else % If its coordinate if larger than 0.5, a negative variation is applied
                table_outputs(n,1+i) = studied_function([points(n,1:i-1) points(n,i)-variation points(n,i+1:nfac)]);
            end
            table_ee(n,i) = abs(table_outputs(n,1+i)-table_outputs(n,1))/variation; % Elementary effect of the i-th factor at the n-th sampled point.
        else % If the factor is over the limit
            table_ee(n,i) = -1; % Its elementary effect is not calculated
        end
    end
end
%% IDENTIFICATION OF THE FACTORS OVER THE LIMIT
if initialization % If the calculations will be foccused on the factors under the limit at the next step
    % Identification of the maxima of the elementary effects of each factor
    resultant = max(table_ee,[],1);
    [resultant_sort, order] = sort(resultant); % Ordering the maxima of the elementary effects of each factor.
    resultant_sort_shift = [resultant_sort(1) resultant_sort(1:nfac-1)]; % New table where the first elementary effect is repeated twice.
    delta = resultant_sort-resultant_sort_shift; % Variation between the successive maxima.
    gagnant = find(delta>=(max(delta)*options.limit_factor), 1 ); % Limit = First variation larger than or equal to the largest variation between two successive maxima
        
    % Is there any missing elementary effect for the factors under the
    % limit ? (= a factor under the limit was formerly over the limit)
    if ~isempty(factors_over) % If there were factors over the limit after the last step
        factors_problem = []; % Factors with missing elementary effects
        for i=1:length(factors_over) % For each of the factors which were over the limit after the last step
            indice_autre_facteur = find(order==factors_over(i)); % Index of the factor in the new ordered values.
            if indice_autre_facteur<gagnant % If the factor is before the limit
                factors_problem = [factors_problem factors_over(i)]; % The list of factors with missing elementary effects is added.
            end
        end
        if ~isempty(factors_problem) % If one of the factors which were over the limit is now under the limit
            
            verification_factors = 0; % Stopping variable.
            i = 1; % From the first step.
            while verification_factors==0 && i<=size(table_ee,1) % While the stopping variable has not been enabled, and while there are still calculation steps to check.
                verification_factors = sum(table_ee(i,factors_problem)==-1); % If there is no elementary effect for one or more of the factors at the i-th step, verification_factors becomes different from zero
                if ~verification_factors % If the stopping variable is still not enabled, the next step is considered
                    i=i+1;
                end
            end
            if verification_factors % If the stopping variable has been enabled, it is necessary to return to the step with a missing elementary effect
                n=i;
            else % Or the algorithms continues normally
                n=n+1;
            end
        else % If none of the factors which were over the limit moved under the limit
            n=n+1;
        end
    else % If there were no factor over the limit after the last step
        n=n+1;
    end
    gagnants = gagnant:nfac; % Ordered indexes of the factors over the limit
    factors_over = order(gagnants); % Real indexes of the factors over the limit
else % If the elementary effects will still be calculated for all the factors at the next step
    n=n+1;
    factors_over=[];
end

end