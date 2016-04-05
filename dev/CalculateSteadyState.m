%%% a function for calculating the steady state density matrix of a given 
%%% set of operators

function rho_ss = CalculateSteadyState( input, init_rho, solution )

    %% machine tolerance error for calculating the zeros of the svd 
    % decomposition
    machine_tol = 1e-10;
    converged = false;

    %% get useful input values
    L = input.L;
    M = input.M;
    solution.rho = init_rho;

    %% create the linear superoperators
    for i = 1:numel(L)
        if isequal(L{i}, input.isLinear)
            Lmat0 = CreateSuperoperatorMatrix( L{i}, input, solution );
        end
    end
    
    %% check whether a steady state solution exists, if not, throw an error
    [steadystateexists, errortype] = CheckForSteadyState(Lmat0);
    if ~steadystateexists
        exception = MException('CalculateSteadyState:InputInvalid', ...
            strcat('the system described produces a non-relaxing dynamical map: ', errortype));
        throw(exception)
    end
        
    %% create the free evolution superoperator matrix and calculate its
    % eigenvalues and eigenvectors
    if ~isfield(solution, 'Uinv')
        Lmat0 = CreateSuperoperatorMatrix(@L0, input, solution);
        [U, D] = eig(Lmat0);
        Uinv = eye(M^2) / U;
        
        solution.U = U;
        solution.Uinv = Uinv;
        solution.D = sparse(D);
    end
    
    %% perform the steady state calculation iteration
    for i = 1:input.Niter

        %% create the new superoperators matrices
        Liter = Lmat0;
        for j = 1:numel(L)
            if ~isequal(L{j}, input.isLinear)
                Liter = Liter + CreateSuperoperatorMatrix(L{j}, input, solution);
            end
        end
        
        %% calculate the SVD of Liter
        [U, S, V] = svd(Liter);
        
        %% the singular value is the solution to Liter * new_rho = 0
        k = find(diag(S) < machine_tol, 1);
        new_rho = V(:, k);
        
        %% reshape and normalise
        new_rho = reshape(new_rho, [M M]);
        new_rho = new_rho / trace(new_rho);
        
        %% check whether the solution has converged
        error = max(max(abs(new_rho - solution.rho)));
        if error < input.SSError
            converged = true;
            break;
        end
        
        %% set the solution density matrix to be this new one
        solution.rho = new_rho;
        
    end
    
    %% check to make sure the solution is converged
    if ~converged
        exception = MException('CalculateSteadyState:NotConverged', ...
            'convergence not acheived');
        throw(exception)
    end
    
	%% return the calculated value    
    solution.rho = reshape(solution.rho, [M M]);
    rho_ss = solution.rho;

end

