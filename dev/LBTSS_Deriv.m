%%% a function to create the Jacobian of the superoperator matrix for the 
%%% steady state Born term operator

function result = LBTSS_Deriv( input, mat, solution )

    %% extract the interactions and other data from the input
    M = input.M;
    interactions = input.interactions;
    partitionIndex = input.partitionIndex;
    rho = solution.rho;
    
    %% create the free evolution superoperator matrix and calculate its
    % eigenvalues and eigenvectors if not already done
    if ~isfield(solution, 'Uinv')
        Lmat0 = CreateSuperoperatorMatrix(@L0, input, solution);
        [U, D] = eig(full(Lmat0));
        Uinv = eye(M^2) / U;
       
        solution.U{partitionIndex} = U;
        solution.Uinv{partitionIndex} = Uinv;
        solution.D{partitionIndex} = sparse(D);
    end
    
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
    LBTSS_deriv = zeros(M^2, M^2);
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
        
            % calculate fluctuations using the previous iteration
            % density matrix in the traces
            dAk = Ak.Operator - trace(Ak.Operator * rho{partitionIndex}) * speye(size(Ak.Operator));
            dBk = Bk.Operator - trace(Bk.Operator * rho{Bk.Index}) * speye(size(Bk.Operator));
            
            %% cycle through the projection matrices
            for count = 1:M^2
                % create sampling matrices
                Mcount = zeros(M^2, 1);
                Mcount(count) = 1.0;
                Mcount = reshape(Mcount, [M M]);
                
                % calculate the fluctuation derivative
                dAk_deriv = -1.0 * trace(Ak.Operator * Mcount) * speye(size(Ak.Operator));
                dBk_deriv = -1.0 * trace(Bk.Operator * Mcount) * speye(size(Bk.Operator));
                
                rhs = zeros(M, M);
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
                        
                        % and their derivatives
                        d_deriv = trace(Bl.Operator * reshape(Mb * reshape(dBk_deriv * rho{Bk.Index} + dBk * Mcount, [M^2 1]), [M M]));
                        s_deriv = trace(Bl.Operator * reshape(Mb * reshape(rho{Bk.Index} * dBk_deriv + Mcount * dBk, [M^2 1]), [M M]));

                        % calculate the commutator between the steady state 
                        % system density matrix and the fluctuations dAk
                        tmprhs = (dAk_deriv * rho{partitionIndex} + dAk * Mcount) * d + dAk * rho{partitionIndex} * d_deriv ...
                                - (rho{partitionIndex} * dAk_deriv + Mcount * dAk) * s - rho{partitionIndex} * dAk * s_deriv;

                        % operate on this commutator with Ma and weight
                        tmprhs = -1.0 * Dinv(a, b) * Ma * reshape(tmprhs, [M^2 1]);

                        % add this to the right hand side
                        rhs = rhs + reshape(tmprhs, [M M]);
                    end
                end

                %% and evaluate the commutator for this set of interactions
                LBTSS_deriv(:, count) = LBTSS_deriv(:, count) ...
                                        - J2 * reshape(Al.Operator * rhs - rhs * Al.Operator, [M^2 1]);

            end
        end
    end
    
    result = LBTSS_deriv;

end

