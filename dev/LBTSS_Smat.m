%%% a function to create the superoperator matrix for the steady state Born
%%% term operator

function result = LBTSS_Smat( input, solution )

    %% extract the interactions and other data from the input
    M = input.M;
    interactions = input.interactions;
    partitionIndex = input.partitionIndex;
    rho = solution.rho;
    
    %% create prototypes for the superoperator matrix product functions
    % matprodleft and matprodright
    lintmp = reshape(1:M^2, [M M]);
    protoleft = kron(speye(M), lintmp);
    protoright = kron(lintmp', speye(M));
    
    %% precalculate the projector matrices Ma and Mb
    Mmat = cell(M^2, 1);
    for a = 1:M^2
        Mmat{a} = solution.U{partitionIndex} * sparse(a, a, 1.0, M^2, M^2) * solution.Uinv{partitionIndex};
    end
    %% and reciprocals of the eigenvalues
    Dinv = zeros(M^2);
    for a = 1:M^2
        for b = 1:M^2
            Dinv(a, b) = 1.0 / (solution.D{partitionIndex}(a, a) + solution.D{partitionIndex}(b, b));
        end
    end
    
    %% calculate the Born term commutator
    LBTSS = sparse(M^2, M^2);
    %% cycle through each of the terms in the interaction Hamiltonian
    for k = 1:numel(interactions)
        for l = 1:numel(interactions)
            
            if ~interactions{k}.Correlations(l)
                continue
            end
            
            %% get the relevant system operators
            J2 = interactions{k}.interactionStrength * interactions{l}.interactionStrength;
            Ak = interactions{k}.A;
            Bk = interactions{k}.B;
            Al = interactions{l}.A;
            Bl = interactions{l}.B;
            % check whether there is information about redundancies in the
            % interactions which allows us to not recalculate some terms in
            % the Born term sum
            if isfield(interactions{k}, 'coordination') && isfield(interactions{l}, 'coordination')
                if interactions{k}.coordination ~= interactions{l}.coordination
                    exception = MException('LBT:InputInvalid', ...
                        'coordination number of correlated Born terms must be the same');
                    throw(exception)
                end
                
                coordination = interactions{k}.coordination;
            else
                coordination = 1;
            end
        
            % calculate fluctuations using the previous iteration
            % density matrix in the traces
            dAk = Ak.Operator - trace(Ak.Operator * rho{partitionIndex}) * speye(size(Ak.Operator));
            dBk = Bk.Operator - trace(Bk.Operator * rho{Bk.Index}) * speye(size(Bk.Operator));
            
            % precalculate the left and right matrix multiplications of dAk
            dAk_left = matprod(dAk, protoleft);
            dAk_right = matprod(dAk, protoright);
            
            rhs = sparse(M^2, M^2);
            %% calculate the superoperator projector sum
            for a = 1:M^2
                % calculate the projectors for the free evolution
                % superoperator matrix
                Ma = Mmat{a};
                for b = 1:M^2
                    % skip if a or b are zero
                    if abs(solution.D{partitionIndex}(a, a)) < 1e-10 || abs(solution.D{partitionIndex}(b, b)) < 1e-10
                        continue;
                    end
                    % calculate Mb free evolution
                    Mb = Mmat{b};
                    
                    % calculate the memory kernel integrand weights
                    d = trace(Bl.Operator * reshape(Mb * reshape(dBk * rho{Bk.Index}, [M^2 1]), [M M]));
                    s = trace(Bl.Operator * reshape(Mb * reshape(rho{Bk.Index} * dBk, [M^2 1]), [M M]));
                    
                    % calculate the commutator between the steady state 
                    % system density matrix and the fluctuations dAk
                    tmprhs = dAk_left * d - dAk_right * s;
                    
                    % operate on this commutator with Ma and weight
                    tmprhs = -1.0 * Dinv(a, b) * Ma * tmprhs;
                    
                    % add this to the right hand side
                    rhs = rhs + tmprhs;
                end
            end
            
            %% and evaluate the commutator for this set of interactions
            LBTSS = LBTSS - J2 * coordination * ...
                            (matprod(Al.Operator, protoleft) * rhs - matprod(Al.Operator, protoright) * rhs);
            
        end
    end
    
    result = LBTSS;

end

