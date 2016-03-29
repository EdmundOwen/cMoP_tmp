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

    %% perform density matrix iteration
    for i = 1:input.Nt
        
        if input.memoryNeeded
            % save results from the old time step 
            result = SaveMemory(result, input, i);
            result.hist{i}.rho = result.rho;
            
            % iterate the memory functions to the new time
            result = IterateMemory(result, input);
        end
        
        % iterate the density matrix
        result.rho = PerformTimeStep(result.rho, result, L, input);
        
    end
    
end