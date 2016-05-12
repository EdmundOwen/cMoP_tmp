%%% A function to save the nearest-neighbout correlations on each site.  
%%% Only called if needed

function result = SaveCorrelation( result, input, i )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    
    % recover the Hilbert space structure
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    % create the annihilation operator
    a = annihilation(onsitedim);
    
    % calculate the nearest-neighbour correlations for each site
    g2 = zeros(clustersize-1, 1);
    for j = 1:clustersize-1
        g2op = InsertOperator(kron(a, a'), j, onsitedim, clustersize-1);
        g2(j) = trace(g2op * rho);
    end

    % save the data
    result.hist{i}.g2 = g2;
    
end