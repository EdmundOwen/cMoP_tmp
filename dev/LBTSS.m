%%% steady state Born term superoperator

function [result, solution] = LBTSS( input, mat, solution )

    %% extract the interactions and other data from the input
    M = input.M;
    interactions = input.interactions;
    rho = solution.rho;
        
    %% create the free evolution superoperator matrix and calculate its
    % eigenvalues and eigenvectors if not already done
    if ~isfield(solution, 'Uinv')
        Lmat0 = CreateSuperoperatorMatrix(@L0, input, solution);
        [U, D] = eig(full(Lmat0));
        Uinv = eye(M^2) / U;
       
        solution.U = U;
        solution.Uinv = Uinv;
        solution.D = sparse(D);
    end
        
    %% calculate the Born term commutator
    LBTSS = zeros(M, M);
    %% cycle through each of the terms in the interaction Hamiltonian
    for k = 1:numel(interactions)
        for l = 1:numel(interactions)
            
            if ~(((k == 1 || k == 2) && (l == 1 || l == 2)) || ((k == 3 || k == 4) && (l == 3 || l == 4)))
                continue
            end
            
            %% get the relevant system operators
            J2 = interactions{k}{1} * interactions{l}{1};
            Al = interactions{l}{2};
            Bl = interactions{l}{3};
            Ak = interactions{k}{2};
            Bk = interactions{k}{3};            
        
            % calculate fluctuations using the previous iteration
            % density matrix in the traces
            dAk = Ak - trace(Ak * rho) * speye(size(Ak));
            dBk = Bk - trace(Bk * rho) * speye(size(Bk));
            
            rhs = zeros(M, M);
            %% calculate the superoperator projector sum
            for a = 1:M^2
                % calculate the projectors for the free evolution
                % superoperator matrix
                Ma = solution.U * sparse(a, a, 1.0, M^2, M^2) * solution.Uinv;
                for b = 1:M^2
                    % skip if a or b are zero
                    if abs(solution.D(a, a)) < 1e-10 || abs(solution.D(b, b)) < 1e-10
                        continue;
                    end
                    % calculate Mb free evolution
                    Mb = solution.U * sparse(b, b, 1.0, M^2, M^2) * solution.Uinv;
                    
                    % calculate the memory kernel integrand weights
                    d = trace(Bl * reshape(Mb * reshape(dBk * rho, [M^2 1]), [M M]));
                    s = trace(Bl * reshape(Mb * reshape(rho * dBk, [M^2 1]), [M M]));
                    
                    % calculate the commutator between the steady state 
                    % system density matrix and the fluctuations dAk
                    tmprhs = dAk * mat * d - mat * dAk * s;
                    
                    % operate on this commutator with Ma and weight
                    tmprhs = -1.0 / (solution.D(a, a) + solution.D(b, b)) ...
                                       * Ma * reshape(tmprhs, [M^2 1]);
                    
                    % reshape and add this to the right hand side
                    rhs = rhs + reshape(tmprhs, [M M]);
                end
            end
            
            %% and evaluate the commutator for this set of interactions
            LBTSS = LBTSS - J2 * (Al * rhs - rhs * Al);
            
        end
    end
    
    result = LBTSS;

end

