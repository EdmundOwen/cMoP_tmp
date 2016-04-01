%%% a heun iteration scheme (predictor-corrector method)
function [mat, solution] = Heun(mat, solution, L, input, dt)
 
    rhs1 = zeros(size(mat));
    rhs2 = zeros(size(mat));

    %% calculate the rhs of the Loiuville equation
    for j = 1:numel(L)
        rhs1 = rhs1 + L{j}( input, mat, solution );
    end

    %% calculate the euler iteration to make the initial guess
    % (this is the predictor step)
    mat_tilde = mat + dt * rhs1;
    
    %% iterate the memory functions (if needed)
    if input.memoryNeeded && ~input.isMemory
        solution = IterateMemory(solution, input, dt);
    end
    
    %% recalculate the rhs of the Louiville equation using the predicted values
    for j = 1:numel(L)
        rhs2 = rhs2 + L{j}( input, mat_tilde, solution );
    end
    
    %% recalculate the euler iteration assuming that the rhs is the average
    % value (this is the corrector step)
    mat = mat + 0.5 * dt * (rhs1 + rhs2);
    
end