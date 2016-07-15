%%% steady state Born term superoperator

function result = LBTSS( input, mat, solution )

    %% extract the interactions and other data from the input
    M = input.M;
    partitionIndex = input.partitionIndex;
        
    %% create the free evolution superoperator matrix and calculate its
    % eigenvalues and eigenvectors if not already done
    if ~isfield(solution, 'Uinv')
        Lmat0 = CreateSuperoperatorMatrix(@L0, input, solution);
        [U, D] = eig(full(Lmat0));
        Uinv = eye(M^2) / U;
       
        solution.U{partitionIndex} = U;
        solution.Uinv{partitionIndex} = Uinv;
        solution.D{partitionIndex} = sparse(D);
    end
    
    %% get the superoperator matrix
    Lmat = LBTSS_Smat(input, solution);
    
    %% reshape the input matrix and output the corresponding density matrix
    result = Lmat * reshape(mat, [M^2 1]);
    result = reshape(result, [M M]);
    
end

