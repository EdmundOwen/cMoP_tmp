%%% A function to iterate the memory functions using a simple Euler
%%% iteration scheme

function result = IterateMemory( result, input, dt )

    % create flag to tell the iteration process that this is a memory
    % matrix
    input.isMemory = true;

    % the memory evolution is free
    L = { @L0 };

    % iterate for each previous time step
    for j = 1:numel(result.hist)
        for k = 1:numel(input.interactions)
            result.hist{j}.c{k} = PerformTimeStep(result.hist{j}.c{k}, result, L, input, dt);
            result.hist{j}.dkern{k} = PerformTimeStep(result.hist{j}.dkern{k}, result, L, input, dt);
            result.hist{j}.r{k} = PerformTimeStep(result.hist{j}.r{k}, result, L, input, dt);
            result.hist{j}.skern{k} = PerformTimeStep(result.hist{j}.skern{k}, result, L, input, dt);
        end
    end
    
    % reset the memory flag
    input.isMemory = false;
    
end

