%%% a heun iteration scheme (predictor-corrector method)
function [mat, solution] = Heun(mat, solution, L, input, dt)
 
    %% allocate matrices
    k1 = cell(1, input.noPartitions);
    k2 = cell(1, input.noPartitions);
    mat_tilde = cell(1, input.noPartitions);

    %% Calculate values at the initial time point
    for i = 1:input.noPartitions
        % create k matrices
        k1{i} = zeros(size(mat{i}));
        k2{i} = zeros(size(mat{i}));

        %% calculate the rhs of the Loiuville equation
        for j = 1:numel(L)
            k1{i} = k1{i} + L{j}( input.subinput{i}, mat{i}, solution );
        end

        %% calculate the euler iteration to make the initial guess
        % (this is the predictor step)
        mat_tilde{i} = mat{i} + dt * k1{i};
    end
    
    %% iterate the memory functions (if needed)
    for i = 1:input.noPartitions
        if input.memoryNeeded && ~input.isMemory
            solution = IterateMemory(solution, input.subinput{i}, dt, 'heun');
        end
    end
    
    %% update the density matrix in solution
    if ~input.isMemory
        for i = 1:input.noPartitions
            solution.rho{i} = solution.rho{i} + dt * k1{i};
        end
    end
    
    %% Calculate values at the initial time point
    for i = 1:input.noPartitions
        %% recalculate the rhs of the Louiville equation using the predicted values
        for j = 1:numel(L)
            k2{i} = k2{i} + L{j}( input.subinput{i}, mat_tilde{i}, solution );
        end

        %% recalculate the euler iteration assuming that the rhs is the average
        % value (this is the corrector step)
        mat{i} = mat{i} + 0.5 * dt * (k1{i} + k2{i});
    end
    
end