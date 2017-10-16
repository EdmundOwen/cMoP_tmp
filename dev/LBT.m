%%% Born term superoperator

function result = LBT( input, mat, solution )

    %% assert that rho is a struct. otherwise there are no memory functions
    if ~isa(solution, 'struct')
        exception = MException('LBT:InputInvalid', ...
            'input solution must be a struct');
        throw(exception)
    end
    
    %% extract the interactions and other data from the input
    M = input.M;
    dt = input.dt;
    interactions = input.interactions;
    partitionIndex = input.partitionIndex;
    
    %% calculate the Born term commutator
    LBT = zeros(M, M);
    
    %% cycle through each of the terms in the interaction Hamiltonian
    for k = 1:numel(interactions)
        for l = 1:numel(interactions)
            
            if ~interactions{k}.Correlations(l)
                continue
            end
            
            %% get the relevant system operators
            % and calculate the correct interaction constant depending on the type of operators            
            J2 = interactions{k}.interactionStrength * interactions{l}.interactionStrength;
            int_type_k = GetFromInput(interactions{k}, 'interactionType', 'unitary');
            int_type_l = GetFromInput(interactions{l}, 'interactionType', 'unitary');
            if strcmp(int_type_k, 'unitary') && strcmp(int_type_l, 'unitary')
                ;
            elseif (strcmp(int_type_k, 'dissipative') && strcmp(int_type_l, 'unitary')) ...
                        || (strcmp(int_type_k, 'unitary') && strcmp(int_type_l, 'dissipative'))
                J2 = -1i * J2;
            elseif strcmp(int_type_k, 'dissipative') && strcmp(int_type_l, 'dissipative')
                J2 = -J2;
            else
                throw exception
            end
            Al = interactions{l}.A;
            Bl = interactions{l}.B;
            if interactions{k}.coordination ~= interactions{l}.coordination
                exception = MException('LBT:InputInvalid', ...
                    'coordination number of correlated Born terms must be the same');
                throw(exception)
            end
            coordination = interactions{k}.coordination;
        
            rhs = zeros(M, M);
            %% perform rhs integration
            for tprime = 1:numel(solution.hist{partitionIndex})
            
                %% calculate the slice with the contracted d and s term
                if strcmp(int_type_l, 'dissipative_ARB')
                    slice = solution.hist{partitionIndex}{tprime}.c{k} ...
                                            * trace(Bl.Operator' * solution.hist{partitionIndex}{tprime}.dkern{k}) ...
                             - solution.hist{partitionIndex}{tprime}.r{k} ...
                                            * trace(Bl.Operator' * solution.hist{partitionIndex}{tprime}.skern{k});
                else
                    slice = solution.hist{partitionIndex}{tprime}.c{k} ...
                                            * trace(Bl.Operator * solution.hist{partitionIndex}{tprime}.dkern{k}) ...
                             - solution.hist{partitionIndex}{tprime}.r{k} ...
                                            * trace(Bl.Operator * solution.hist{partitionIndex}{tprime}.skern{k});
                end
            
                %% add the slice the integral using the trapezium rule
                rhs = rhs + 0.5 * dt * slice;
                if tprime ~= 1 || tprime ~= numel(solution.hist{partitionIndex})
                    rhs = rhs + 0.5 * dt * slice;
                end
            end
        
            %% and evaluate the commutator for this set of interactions
            if strcmp(int_type_l, 'dissipative_BRA')
                LBT = LBT + J2 * coordination * (Al.Operator' * rhs - rhs * Al.Operator');
            else
                LBT = LBT - J2 * coordination * (Al.Operator * rhs - rhs * Al.Operator);
            end
        
        end
    end
    
    result = LBT;

end

