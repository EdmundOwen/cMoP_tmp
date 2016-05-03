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

    %% create dummy new_rho
    new_rho = result.rho;
    
    %% perform density matrix iteration
    for i = 1:input.Nt
        
        result.rho = new_rho;
        if mod(i, input.NTest) == 0
            for j = 1:numel(input.probelist)
                % save results from the old time step 
                result = input.probelist{j}(result, input, i);
            end
        end
        
        % iterate the density matrix and memory functions
        [new_rho, result] = PerformTimeStep(result.rho, result, L, input, dt);
        
    end
    
    %% and save the result from the final time step
    result.rho = new_rho;
    
end