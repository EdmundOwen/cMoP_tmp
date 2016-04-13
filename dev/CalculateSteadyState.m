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
        [U, D] = eig(full(Lmat0));
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

