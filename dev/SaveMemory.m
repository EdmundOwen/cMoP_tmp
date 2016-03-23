%%% A function to save the memory functions.  Only called if needed

function result = SaveMemory( result, input, i )
        
    for k = 1:numel(input.interactions)
        
        % calculate the fluctuations
        Ak = input.interactions{k}{2}; Bk = input.interactions{k}{3};
        dAk = Ak - trace(Ak * result.rho) * speye(size(Ak));
        dBk = Bk - trace(Bk * result.rho) * speye(size(Bk));
        
        % save the memory functions from this time step
        result.hist{i}.c{k} = dAk * result.rho;
        result.hist{i}.dkern{k} = dBk * result.rho;
        result.hist{i}.r{k} = result.rho * dAk;
        result.hist{i}.skern{k} = result.rho * dBk;
    end

end

