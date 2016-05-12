%%% a function to form the tensor sum of the local Hamiltonians

function H = TensorAddOperator( Hloc, onsitedim, clustersize )

    if clustersize - 1 < 0
        H = 0;
        return
    end

    % initialise the Hamiltonian
    H = kron(Hloc, speye(onsitedim^(clustersize-1)));
    % and add the rest
    for i = 2:clustersize
        H = H + kron(speye(onsitedim^(i - 1)), kron(Hloc, speye(onsitedim^(clustersize - i))));
    end

end

