%%% A function to iterate the memory functions using a simple Euler
%%% iteration scheme

function result = IterateMemory( result, input )

    % the memory evolution is free
    L = { @L0 };

    % iterate for each previous time step
    for j = 1:numel(result.hist)
        for k = 1:numel(input.interactions)
            result.hist{j}.c{k} = PerformTimeStep(result.hist{j}.c{k}, result, L, input);
            result.hist{j}.dkern{k} = PerformTimeStep(result.hist{j}.dkern{k}, result, L, input);
            result.hist{j}.r{k} = PerformTimeStep(result.hist{j}.r{k}, result, L, input);
            result.hist{j}.skern{k} = PerformTimeStep(result.hist{j}.skern{k}, result, L, input);
        end
    end
end

