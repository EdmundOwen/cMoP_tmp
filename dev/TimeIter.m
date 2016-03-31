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
        end
        
        % iterate the density matrix and memory functions
        [new_rho, result] = PerformTimeStep(result.rho, result, L, input, dt);
        result.rho = new_rho;
        
    end
    
end