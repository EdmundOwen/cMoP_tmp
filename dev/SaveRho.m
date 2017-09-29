%%% A function to save the density matrox.  Only called if needed

function solution = SaveRho( solution, input, timestep )
        
    for i = 1:input.noPartitions
        % get the density matrix
        subinput = input.subinput{i};
        partitionIndex = subinput.partitionIndex;
        rho = solution.rho{partitionIndex};
    
        % save this to the history
        solution.hist{partitionIndex}{timestep}.rho = rho;
    end

end

