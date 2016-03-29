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
    
    %% calculate the Born term commutator
    LBT = zeros(M, M);
    
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
        
            rhs = zeros(M, M);
            %% perform rhs integration
            for tprime = 1:numel(solution.hist)
            
                %% calculate the slice with the contracted d and s term
                slice = solution.hist{tprime}.c{k} * trace(Bl * solution.hist{tprime}.dkern{k}) ...
                            - solution.hist{tprime}.r{k} * trace(Bl * solution.hist{tprime}.skern{k});
            
                %% add the slice the integral
                rhs = rhs + dt * slice;
            end
        
            %% and evaluate the commutator for this set of interactions
            LBT = LBT - J2 * (Al * rhs - rhs * Al);
        
        end
    end
    
    result = LBT;

end

