%%% A function to save the density on each site.  Only called if needed

function solution = SaveDensity( solution, input, timestep )

    for i = 1:input.noPartitions
        % get the density matrix 
        subinput = input.subinput{i};
        partitionIndex = subinput.partitionIndex;
        rho = solution.rho{partitionIndex};
    
        % recover the Hilbert space structure
        onsitedim = subinput.onsitedim;
        clustersize = subinput.clustersize;
    
        % create the annihilation operator
        a = annihilation(onsitedim);
    
        % calculate the density on each site
        n = zeros(clustersize, 1);
        for j = 1:clustersize
            nop = InsertOperator(a' * a, j, onsitedim, clustersize);
            n(j) = abs(trace(nop * rho));
        end

        % save the density data
        solution.hist{partitionIndex}{timestep}.n = n;
    end
    
end