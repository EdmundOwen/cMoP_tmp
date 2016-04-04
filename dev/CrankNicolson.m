%%% a crank-nicolson scheme

function [mat, solution] = CrankNicolson(mat, solution, L, input, dt)

    % ensure that M is correct
    input.M = input.onsitedim^input.clustersize;
    
    %% create the right-hand side matrix of the crank nicolson method
    solution = updateLmatrices( L, solution, input );
    rhs = eye(input.M^2);
    for i = 1:numel(L)
        rhs = rhs + 0.5 * dt * solution.Lmat{i};
    end
    
    %% iterate the memory functions (if needed)
    if input.memoryNeeded && ~input.isMemory
        solution = IterateMemory(solution, input, dt);
    end
    
    %% create the right-hand side matrix of the crank nicolson method
    solution = updateLmatrices( L, solution, input );
    lhs = eye(input.M^2);
    for i = 1:numel(L)
        lhs = lhs - 0.5 * dt * solution.Lmat{i};
    end
    
    %% solve and reshape
    b = rhs * reshape(mat, [input.M^2 1]);
    a = lhs \ b;
    
    mat = reshape(a, [input.M input.M]);
    
end

%%% update L matrices
function solution = updateLmatrices( L, solution, input )
    
    %% create the cell array for the L matrices
    if ~isfield(solution, 'Lmat') && ~input.isMemory
        solution.Lmat = cell([numel(L) 1]);
    end
    
    %% cycle through superoperators
    % if the matrix doesn't exist and superoperator doesn't need memory terms;
    % or if the term uses memory, save it to solution
    for i = 1:numel(L)
        for j = 1:numel(input.termsNeedingMemory)
            if isempty(solution.Lmat{i}) && ~isequal(L{i}, input.termsNeedingMemory{j}) || ...
                                                        isequal(input.termsNeedingMemory{j}, L{i})
                
                solution.Lmat{i} = CreateSuperoperatorMatrix(L{i}, input, solution);
                
            end
        end
    end
    
end
