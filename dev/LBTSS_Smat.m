%%% a function to create the superoperator matrix for the steady state Born
%%% term operator

function result = LBTSS_Smat( input, solution )

    %% extract the interactions and other data from the input
    M = input.M;
    interactions = input.interactions;
    rho = solution.rho;
    
    %% create prototypes for the superoperator matrix product functions
    % matprodleft and matprodright
    lintmp = reshape(1:M^2, [M M]);
    protoleft = kron(speye(M), lintmp);
    protoright = kron(lintmp', speye(M));
    
    %% calculate the Born term commutator
    LBTSS = sparse(M^2, M^2);
    %% cycle through each of the terms in the interaction Hamiltonian
    for k = 1:numel(interactions)
        for l = 1:numel(interactions)
            
            if ~interactions{k}{4}(l)
                continue
            end
            
            %% get the relevant system operators
            J2 = interactions{k}{1} * interactions{l}{1};
            Ak = interactions{k}{2};
            Bk = interactions{k}{3};
            Al = interactions{l}{2};
            Bl = interactions{l}{3};
        
            % calculate fluctuations using the previous iteration
            % density matrix in the traces
            dAk = Ak - trace(Ak * rho) * speye(size(Ak));
            dBk = Bk - trace(Bk * rho) * speye(size(Bk));
            
            rhs = sparse(M^2, M^2);
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
                    tmprhs = matprod(dAk, protoleft) * d - matprod(dAk, protoright) * s;
                    
                    % operate on this commutator with Ma and weight
                    tmprhs = -1.0 / (solution.D(a, a) + solution.D(b, b)) ...
                                       * Ma * tmprhs;
                    
                    % add this to the right hand side
                    rhs = rhs + tmprhs;
                end
            end
            
            %% and evaluate the commutator for this set of interactions
            LBTSS = LBTSS - J2 * (matprod(Al, protoleft) * rhs - matprod(Al, protoright) * rhs);
            
        end
    end
    
    result = LBTSS;

end

