%%% A function to save the density on each site.  Only called if needed

function result = SaveDensity( result, input, timestep )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    partitionIndex = input.partitionIndex;
    
    % recover the Hilbert space structure
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    % create the annihilation operator
    a = annihilation(onsitedim);
    
    % calculate the magnetisation on each site
    sigmax = zeros(clustersize, 1);
    sigmay = zeros(clustersize, 1);
    sigmaz = zeros(clustersize, 1);
    for j = 1:clustersize
        sigmax_op = InsertOperator(a + a', j, onsitedim, clustersize);
        sigmay_op = InsertOperator(1i * (a - a'), j, onsitedim, clustersize);
        sigmaz_op = InsertOperator(a' * a - a * a', j, onsitedim, clustersize);
        sigmax(j) = real(trace(sigmax_op * rho{partitionIndex}));
        sigmay(j) = real(trace(sigmay_op * rho{partitionIndex}));
        sigmaz(j) = real(trace(sigmaz_op * rho{partitionIndex}));
    end

    % save the density data
    result.hist{partitionIndex}{timestep}.sigmax = sigmax;
    result.hist{partitionIndex}{timestep}.sigmay = sigmay;
    result.hist{partitionIndex}{timestep}.sigmaz = sigmaz;
end