%%% A function to iterate the memory functions using a simple Euler
%%% iteration scheme

function result = IterateMemory( result, input )

    dt = input.dt;

    % iterate for each previous time step
    for j = 1:numel(result.hist)
        for k = 1:numel(input.interactions)
            result.hist{j}.c{k} = result.hist{j}.c{k} + dt * L0( input, result.hist{j}.c{k} );
            result.hist{j}.dkern{k} = result.hist{j}.dkern{k} + dt * L0( input, result.hist{j}.dkern{k} );
            result.hist{j}.r{k} = result.hist{j}.r{k} + dt * L0( input, result.hist{j}.r{k} );
            result.hist{j}.skern{k} = result.hist{j}.skern{k} + dt * L0( input, result.hist{j}.skern{k} );
        end
    end

end

