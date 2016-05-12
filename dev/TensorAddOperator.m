%%% a function to form the tensor sum of the local Hamiltonians

function H = TensorAddOperator( Hloc, onsitedim, clustersize )

    if clustersize - 1 < 0
        H = 0;
        return
    end

    % initialise the Hamiltonian
    H = InsertOperator(Hloc, 1, onsitedim, clustersize);
    % and add the rest
    for i = 2:clustersize
        H = H + InsertOperator(Hloc, i, onsitedim, clustersize);
    end

end

