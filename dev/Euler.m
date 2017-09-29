%%% a simple euler iteration scheme

function [mat, solution] = Euler(mat, solution, L, input, dt)

    for i = 1:input.noPartitions
        
        rhs = zeros(size(mat{i}));
        %% calculate rhs of the Louiville equation
        for j = 1:numel(L)
            rhs = rhs + L{j}( input.subinput{i}, mat{i}, solution );
        end
  
        %% iterate the matrix
        mat{i} = mat{i} + dt * rhs;
    end
  
    %% iterate the memory functions (if needed)
    for i = 1:input.noPartitions
        if input.memoryNeeded && ~input.isMemory
            solution = IterateMemory(solution, input.subinput{i}, dt, 'euler');
        end
    end
    
end