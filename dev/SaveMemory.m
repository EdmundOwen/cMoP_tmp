%%% A function to save the memory functions.  Only called if needed

function result = SaveMemory( result, input, timestep )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    partitionIndex = input.partitionIndex;
    
    for k = 1:numel(input.interactions)
        % calculate the fluctuations
        Ak = input.interactions{k}.A; 
        Bk = input.interactions{k}.B;
        
        dAk = Ak.Operator - trace(Ak.Operator * rho{partitionIndex}) * speye(size(Ak.Operator));
        dBk = Bk.Operator - trace(Bk.Operator * rho{Bk.Index}) * speye(size(Bk.Operator));
        
        % save the memory functions from this time step
        result.hist{partitionIndex}{timestep}.c{k} = dAk * rho{partitionIndex};
        result.hist{partitionIndex}{timestep}.dkern{k} = dBk * rho{Bk.Index};
        result.hist{partitionIndex}{timestep}.r{k} = rho{partitionIndex} * dAk;
        result.hist{partitionIndex}{timestep}.skern{k} = rho{Bk.Index} * dBk;
    end

end

