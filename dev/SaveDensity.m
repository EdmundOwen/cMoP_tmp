%%% A function to save the density on each site.  Only called if needed

function result = SaveDensity( result, input, i )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    
    % recover the Hilbert space structure
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    % create the annihilation operator
    a = annihilation(onsitedim);
    
    % calculate the density on each site
    n = zeros(clustersize, 1);
    for j = 1:clustersize
        nop = InsertOperator(a' * a, j, onsitedim, clustersize);
        n(j) = abs(trace(nop * rho));
    end

    % save the density data
    result.hist{i}.n = n;
end