%%% A function to save the memory functions.  Only called if needed

function result = SaveMemory( result, input, i )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;

    for k = 1:numel(input.interactions)
        % calculate the fluctuations
        Ak = input.interactions{k}{2}; Bk = input.interactions{k}{3};
        dAk = Ak - trace(Ak * rho) * speye(size(Ak));
        dBk = Bk - trace(Bk * rho) * speye(size(Bk));
        
        % save the memory functions from this time step
        result.hist{i}.c{k} = dAk * rho;
        result.hist{i}.dkern{k} = dBk * rho;
        result.hist{i}.r{k} = rho * dAk;
        result.hist{i}.skern{k} = rho * dBk;
    end

end

