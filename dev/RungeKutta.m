%%% a basic runga-kutta iteration scheme
function [mat, solution] = RungeKutta(mat, solution, L, input, dt)
    
    %% allocate k matrices
    k1 = cell(1, input.noPartitions);
    k2 = cell(1, input.noPartitions);
    k3 = cell(1, input.noPartitions);
    k4 = cell(1, input.noPartitions);

    %% Calculate values at the initial time point
    for i = 1:input.noPartitions
        % create k matrices
        k1{i} = zeros(size(mat{i}));
        k2{i} = zeros(size(mat{i}));
        k3{i} = zeros(size(mat{i}));
        k4{i} = zeros(size(mat{i}));

        % calculate rhs of the Louiville equation at the initial time
        for j = 1:numel(L)
            k1{i} = k1{i} + L{j}( input.subinput{i}, mat{i}, solution );
        end
    end

    %% iterate the memory functions half a time step (if needed)
    for i = 1:input.noPartitions
        if input.memoryNeeded && ~input.isMemory
            solution = IterateMemory(solution, input.subinput{i}, 0.5 * dt, 'runge-kutta');
        end
    end
    
    %% Calculate values at the middle time point
    for i = 1:input.noPartitions
        % calculate the rhs at the midpoint using k1
        for j = 1:numel(L)
            k2{i} = k2{i} + L{j}( input.subinput{i}, mat{i} + 0.5 * dt * k1{i}, solution );
        end
        % and using k2
        for j = 1:numel(L)
            k3{i} = k3{i} + L{j}( input.subinput{i}, mat{i} + 0.5 * dt * k2{i}, solution );
        end
    end

    %% iterate the memory functions for the rest of the time step (if needed)
    for i = 1:input.noPartitions
        if input.memoryNeeded && ~input.isMemory
            solution = IterateMemory(solution, input.subinput{i}, 0.5 * dt, 'runge-kutta');
        end
    end
    
    %% Calculate values at the final time point and iterate matrices
    for i = 1:input.noPartitions
        % calculate the rhs at the end of the interval
        for j = 1:numel(L)
            k4{i} = k4{i} + L{j}( input.subinput{i}, mat{i} + dt * k3{i}, solution );
        end

        % iterate the matrix using the runga-kutta values
        mat{i} = mat{i} + (dt / 6.0) * (k1{i} + 2.0 * k2{i} + 2.0 * k3{i} + k4{i});
    end

end

