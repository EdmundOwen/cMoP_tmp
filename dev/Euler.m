%%% a simple euler iteration scheme

function [mat, solution] = Euler(mat, solution, L, input, dt)

    rhs = zeros(size(mat));
    %% calculate rhs of the Louiville equation
    for j = 1:numel(L)
        rhs = rhs + L{j}( input, mat, solution );
    end
  
    %% iterate the memory functions (if needed)
    if input.memoryNeeded && ~input.isMemory
        solution = IterateMemory(solution, input, dt, 'euler');
    end
  
    %% iterate the matrix
    mat = mat + dt * rhs;
    
end