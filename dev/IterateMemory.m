%%% A function to iterate the memory functions using a simple Euler
%%% iteration scheme

function result = IterateMemory( result, input, dt, method )

    % create flag to tell the iteration process that this is a memory
    % matrix
    input.isMemory = true;
    partitionIndex = input.partitionIndex;

    % the memory evolution is free
    L = { @L0 };

    % iterate for each previous time step
    for j = 1:numel(result.hist{partitionIndex})
        for k = 1:numel(input.interactions)
            result.hist{partitionIndex}{j}.c{k} = PerformTimeStep(result.hist{partitionIndex}{j}.c{k}, result, L, input, dt, method);
            result.hist{partitionIndex}{j}.dkern{k} = PerformTimeStep(result.hist{partitionIndex}{j}.dkern{k}, result, L, input, dt, method);
            result.hist{partitionIndex}{j}.r{k} = PerformTimeStep(result.hist{partitionIndex}{j}.r{k}, result, L, input, dt, method);
            result.hist{partitionIndex}{j}.skern{k} = PerformTimeStep(result.hist{partitionIndex}{j}.skern{k}, result, L, input, dt, method);
        end
    end
    
    % reset the memory flag
    input.isMemory = false;
    
end

