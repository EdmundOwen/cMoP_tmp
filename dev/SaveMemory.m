%%% A function to save the memory functions.  Only called if needed

function result = SaveMemory( result, input )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    
    % calculate the length of the history to find out where the memory
    % should be saved to
    hist_length = numel(result.hist);

    for k = 1:numel(input.interactions)
        % calculate the fluctuations
        Ak = input.interactions{k}{2}; Bk = input.interactions{k}{3};
        dAk = Ak - trace(Ak * rho) * speye(size(Ak));
        dBk = Bk - trace(Bk * rho) * speye(size(Bk));
        
        % save the memory functions from this time step
        result.hist{hist_length + 1}.c{k} = dAk * rho;
        result.hist{hist_length + 1}.dkern{k} = dBk * rho;
        result.hist{hist_length + 1}.r{k} = rho * dAk;
        result.hist{hist_length + 1}.skern{k} = rho * dBk;
    end

end

