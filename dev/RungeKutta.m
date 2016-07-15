%%% a basic runga-kutta iteration scheme
function [mat, solution] = RungeKutta(mat, solution, L, input, dt)
    
    k1 = zeros(size(mat));
    k2 = zeros(size(mat));
    k3 = zeros(size(mat));
    k4 = zeros(size(mat));
    
    %% calculate rhs of the Louiville equation at the initial time
    for j = 1:numel(L)
        k1 = k1 + L{j}( input, mat, solution );
    end
    
    %% iterate the memory functions half a time step (if needed)
    if input.memoryNeeded && ~input.isMemory
        solution = IterateMemory(solution, input, 0.5 * dt, 'runge-kutta');
    end
    
    %% calculate the rhs at the midpoint using k1
    for j = 1:numel(L)
        k2 = k2 + L{j}( input, mat + 0.5 * dt * k1, solution );
    end
    %% and using k2
    for j = 1:numel(L)
        k3 = k3 + L{j}( input, mat + 0.5 * dt * k2, solution );
    end

    %% iterate the memory functions for the rest of the time step (if needed)
    if input.memoryNeeded && ~input.isMemory
        solution = IterateMemory(solution, input, 0.5 * dt, 'runge-kutta');
    end
    
    %% calculate the rhs at the end of the interval
    for j = 1:numel(L)
        k4 = k4 + L{j}( input, mat + dt * k3, solution );
    end
    
    %% iterate the matrix using the runga-kutta values
    mat = mat + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    
end

