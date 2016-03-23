%%% Time iteration

function result = TimeIter(input, rho)

    %% check to see if the input is a struct, and if so, extract the density matrix
    if isa(rho, 'struct')
        result.rho = rho.rho;
    else
        result.rho = rho;
    end
    result.hist = {};
    
    %% Get the list of linear operators on the rhs of the Louiville equation
    L = input.L;
    dt = input.dt;

    %% perform density matrix iteration
    for i = 1:input.Nt
        
        if input.memoryNeeded
            % save results from the old time step 
            result = SaveMemory(result, input, i);
            result.hist{i}.rho = result.rho;
            
            % iterate the memory functions to the new time
            result = IterateMemory(result, input);
        end
        
        rhs = zeros(size(result.rho));
        % calculate rhs of Louiville equaion
        for j = 1:numel(L)
            rhs = rhs + L{j}( input, result );
        end
        
        % iterate the density matrix
        result.rho = result.rho + dt * rhs;
        
    end
    
end