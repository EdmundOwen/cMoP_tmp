%%% A function to save the density on each site.  Only called if needed

function solution = SaveMagnetisation( solution, input, timestep )
    
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

        % calculate the magnetisation on each site
        sigmax = zeros(clustersize, 1);
        sigmay = zeros(clustersize, 1);
        sigmaz = zeros(clustersize, 1);
        for j = 1:clustersize
            sigmax_op = InsertOperator(a + a', j, onsitedim, clustersize);
            sigmay_op = InsertOperator(-1i * (a - a'), j, onsitedim, clustersize);
            sigmaz_op = InsertOperator(a * a' - a' * a, j, onsitedim, clustersize);
            sigmax(j) = real(trace(sigmax_op * rho));
            sigmay(j) = real(trace(sigmay_op * rho));
            sigmaz(j) = real(trace(sigmaz_op * rho));
        end

        % save the magnetisation data
        solution.hist{partitionIndex}{timestep}.sigmax = sigmax;
        solution.hist{partitionIndex}{timestep}.sigmay = sigmay;
        solution.hist{partitionIndex}{timestep}.sigmaz = sigmaz;
    end
    
end