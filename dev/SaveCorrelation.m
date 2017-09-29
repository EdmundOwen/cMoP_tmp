%%% A function to save the nearest-neighbout correlations on each site.  
%%% Only called if needed

function result = SaveCorrelation( result, input, timestep )
    
    if (input.noPartitions ~= 1)
            exception = MException('SaveCorrelation:NotImplementedException', ...
                strcat('noPartitions = ', input.noPartitions, ' is not implemented in SaveCorrelation'));
            throw(exception)
    end

    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    partitionIndex = input.partitionIndex;
    
    % recover the Hilbert space structure
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    % create the annihilation operator
    a = annihilation(onsitedim);
    
    % calculate the nearest-neighbour correlations for each site
    g2 = zeros(clustersize-1, 1);
    for j = 1:clustersize-1
        g2op = InsertOperator(kron(a, a'), j, onsitedim, clustersize-1);
        g2(j) = trace(g2op * rho{partitionIndex});
    end

    % save the data
    result.hist{partitionIndex}{timestep}.g2 = g2;
    
end