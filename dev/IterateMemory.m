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
            % create copy of input and set number of partitions to zero
            % (memory evolution is free so there should be no interactions)
            dummyinput = input;
            dummyinput.noPartitions = 1;
            dummyinput.L = L;
            dummyinput.subinput{1} = input;
            
            result.hist{partitionIndex}{j}.c{k} = PerformTimeStep({result.hist{partitionIndex}{j}.c{k}}, result, L, dummyinput, dt, method);
            result.hist{partitionIndex}{j}.dkern{k} = PerformTimeStep({result.hist{partitionIndex}{j}.dkern{k}}, result, L, dummyinput, dt, method);
            result.hist{partitionIndex}{j}.r{k} = PerformTimeStep({result.hist{partitionIndex}{j}.r{k}}, result, L, dummyinput, dt, method);
            result.hist{partitionIndex}{j}.skern{k} = PerformTimeStep({result.hist{partitionIndex}{j}.skern{k}}, result, L, dummyinput, dt, method);
            
            % remove the memory from the cell
            result.hist{partitionIndex}{j}.c{k} = result.hist{partitionIndex}{j}.c{k}{1};
            result.hist{partitionIndex}{j}.dkern{k} = result.hist{partitionIndex}{j}.dkern{k}{1};
            result.hist{partitionIndex}{j}.r{k} = result.hist{partitionIndex}{j}.r{k}{1};
            result.hist{partitionIndex}{j}.skern{k} = result.hist{partitionIndex}{j}.skern{k}{1};
        end
    end
    
    % reset the memory flag
    input.isMemory = false;
    
end

