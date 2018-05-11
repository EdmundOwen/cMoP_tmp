%%% a function for calculating the steady state density matrix of a given 
%%% set of operators

function rho_ss = CalculateSteadyState( input, init_rho, solution )

    %% machine tolerance error for calculating the zeros of the svd 
    % decomposition
    machine_tol = 1e-10;
    alpha = input.steadyStateMixingRatio;
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
        
            %% solutions to Liter * rho = 0 are overconstrained as we require
            % the norm to be conserved.  we use this as the last line of the
            % matrix to solve
            normtrace = reshape(speye(M), [1 M^2]);
            X = sparse([Liter(1:M^2-1, 1:M^2); ...
                        normtrace]);

            % the right hand side of the matrix equation is all zeros apart
            % for the trace preserving part
            B = sparse(M^2, 1, 1.0);

            %% solve and reshape to get the new density matrix
            v = X\B;
            new_rho = reshape(v, [M M]);
        
            %% check whether the solution has converged
            error = max(max(abs(new_rho - solution.rho{k})));
            if error < input.SSError
                if CheckBool(input, 'verbose', false)
                    fprintf('# iterations to convergence: %i\n', i)
                end
                converged(k) = true;
            else
                if CheckBool(input, 'verbose', false)
                    fprintf('Iteration %i:  error = %g\n', i, full(error));
                end
                converged(k) = false;
            end
        
            %% set the solution density matrix to be this new one
            solution.rho{k} = (1.0 - alpha) * solution.rho{k} + alpha * new_rho;
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

