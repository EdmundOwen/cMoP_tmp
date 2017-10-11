%%% a function for calculating the steady state density matrix of a given 
%%% set of operators

function rho_ss = CalculateSteadyStateWithGradient( input, init_rho, solution )

    %% machine tolerance error for calculating the zeros of the svd 
    % decomposition
    machine_tol = 1e-10;
    converged = zeros(1, input.noPartitions);

    %% get useful input values
    L = input.L;
    solution.rho = init_rho;

    Lmat0 = cell(1,input.noPartitions);
    for k = 1:input.noPartitions
        %% create the linear superoperators
        for i = 1:numel(L)
            if strcmp(func2str(L{i}), func2str(input.isLinear))
                Lmat0{k} = CreateSuperoperatorMatrix( L{i}, input.subinput{k}, solution );
            end
        end
    
        %% check whether a steady state solution exists, if not, throw an error
        [steadystateexists, errortype] = CheckForSteadyState(Lmat0{k});
        if ~steadystateexists
            exception = MException('CalculateSteadyState:InputInvalid', ...
                strcat('the system described produces a non-relaxing dynamical map: ', errortype));
            throw(exception)
        end
        
        %% create the free evolution superoperator matrix
        Lmat0{k} = CreateSuperoperatorMatrix(@L0, input.subinput{k}, solution);
    end
    
    %% perform the steady state calculation iteration
    for i = 1:input.Niter

        %% cycle through the constituent systems
        for k = 1:input.noPartitions
            
            M = input.subinput{k}.M;
            %% create the new superoperators matrices
            Liter = Lmat0{k};
            for j = 1:numel(L)
                if ~strcmp(func2str(L{j}), func2str(input.isLinear))
                    Liter = Liter + CreateSuperoperatorMatrix(L{j}, input.subinput{k}, solution);
                end
            end
            
            %% create Jacobian superoperator matrix
            % the free evolution is linear
            LJacob = Lmat0{k};
            if isempty(cell2mat(strfind(cellfun(@func2str, input.subinput{k}.L, 'UniformOutput', 0), 'LMF'))) == 0
                LJacob = LJacob + LMF_Deriv(input.subinput{k}, solution.rho{k}, solution);
            end
            if isempty(cell2mat(strfind(cellfun(@func2str, input.subinput{k}.L, 'UniformOutput', 0), 'LBTSS'))) == 0
                LJacob = LJacob + LBTSS_Deriv(input.subinput{k}, solution.rho{k}, solution);
            end
            
            %% calculate the correction using Newton's method
            %% g' x = -g where g = L * rho
            grho = -1.0 * Liter * reshape(solution.rho{k}, [M^2, 1]);
            
            %% initial solution has trace of unity so x must be traceless
            normtrace = reshape(eye(M), [1 M^2]);
            X = [LJacob(1:M^2-1, 1:M^2); ...
                        normtrace];
            grho(M^2) = 0.0;

            x = reshape(X \ grho, [M M]);
            
            %% for the moment, do binary chop
            alpha1 = input.searchInterval(1);
            alpha2 = input.searchInterval(2);
            searchCount = input.searchCount;
            
            switch input.optimalIncrementSearchMethod
                case 'noSearch'
                    new_rho = solution.rho{k} + x;
                    
                case 'binarychopExact'
                    %% use the exact linear operator to evaluate how much of x to add
                    % for alpha1
                    oldrho = solution.rho{k};
                    solution.rho{k} = oldrho  + alpha1 * x;
                    Lrho = zeros(size(solution.rho{k}));
                    for count = 1:numel(input.L)
                        Lrho = Lrho + L{count}(input.subinput{k}, solution.rho{k}, solution);
                    end
                    feval1 = norm(reshape(Lrho, [M^2 1]));
                    solution.rho{k} = oldrho  + alpha2 * x;
                    % for alpha2
                    Lrho = 0.0 * Lrho;
                    for count = 1:numel(input.L)
                        Lrho = Lrho + L{count}(input.subinput{k}, solution.rho{k}, solution);
                    end                
                    feval2 = norm(reshape(Lrho, [M^2 1]));

                    for count = 1:searchCount
                        solution.rho{k} = oldrho + 0.5 * (alpha1 + alpha2) * x;
                        Lrho = 0.0 * Lrho;
                        for count = 1:numel(input.L)
                            Lrho = Lrho + L{count}(input.subinput{k}, solution.rho{k}, solution);
                        end                
                        fevaltmp = norm(reshape(Lrho, [M^2 1]));

                        if feval2 < feval1
                            alpha1 = 0.5 * (alpha1 + alpha2);
                            feval1 = fevaltmp;
                        else
                            alpha2 = 0.5 * (alpha1 + alpha2);
                            feval2 = fevaltmp;
                        end
                    end
                    % set new_rho to the new value
                    new_rho = oldrho + 0.5 * (alpha1 + alpha2) * x;
                    if input.verbose
                        fprintf('Iteration %i:  alpha = %g\n', i, 0.5 * (alpha1+alpha2));
                    end
                    
                case 'binarychopLinearized'
                    %% use the current linearized version of the superoperator matrix 
                    %% to evaluate how much of x to add (approximate but hopefully faster)
                    feval1 = norm(Liter * reshape(solution.rho{k} + alpha1 * x, [M^2 1]));
                    feval2 = norm(Liter * reshape(solution.rho{k} + alpha2 * x, [M^2 1]));
                    for count = 1:searchCount
                        fevaltmp = norm(Liter * reshape(solution.rho{k} + 0.5 * (alpha1 + alpha2) * x, [M^2 1]));

                        if feval2 < feval1
                            alpha1 = 0.5 * (alpha1 + alpha2);
                            feval1 = fevaltmp;
                        else
                            alpha2 = 0.5 * (alpha1 + alpha2);
                            feval2 = fevaltmp;
                        end
                    end
            
                    % set new_rho to the new value
                    new_rho = solution.rho{k} + 0.5 * (alpha1 + alpha2) * x;
                    
                    if input.verbose
                        fprintf('Iteration %i:  alpha = %g\n', i, 0.5 * (alpha1+alpha2));
                    end
                    
                otherwise
                    exception = MException('CalculateSteadyStateWithGradient:InputInvalid', ...
                        strcat('The method for searching for the optimal amount of x to add for each iteration step (', input.optimalIncrementSearchMethod, ') is not supported: ', errortype));
                    throw(exception)
                    
            end
                    
        
            %% check whether the solution has converged
            error = max(max(abs(new_rho - solution.rho{k})));
            if error < input.SSError
                if CheckBool(input, 'verbose', false)
                    fprintf('# iterations to convergence: %i\n', i)
                end
                converged(k) = true;
            else
%                 toc
                if CheckBool(input, 'verbose', false)
                 %   fprintf('Completed iteration: %i\n', i)
                end
                converged(k) = false;
%                 tic
            end
        
            %% set the solution density matrix to be this new one
            solution.rho{k} = new_rho;
        end
        
        %% check that all the partitions have converged
        if sum(converged) == input.noPartitions
            break;
        end
    end
    
    %% check to make sure the solution is converged
    if sum(converged)~=input.noPartitions
        % throw an exception if these are not being ignored
        if ~CheckBool(input, 'ignoreExceptions', true)
            exception = MException('CalculateSteadyState:NotConverged', ...
                'convergence not acheived');
            throw(exception)
        end
        
        rho_ss = -1;
        return
    end
    
    for k = 1:input.noPartitions
        %% return the calculated value
        M = input.subinput{k}.M;
        solution.rho{k} = reshape(solution.rho{k}, [M M]);
    end
    rho_ss = solution.rho;

end

