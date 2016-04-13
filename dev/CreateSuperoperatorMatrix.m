%%% a function to create superoperator matrices by putting matrices with a 
%%% single one into the superoperator method

function mat = CreateSuperoperatorMatrix( L, input, solution )

    % get the matrix dimension
    M = input.M;

    if isequal(L, @L0)
        mat = L0_Smat(input, solution);
        return
    elseif isequal(L, @LBTSS)
        
        % create the free evolution superoperator matrix and calculate its
        % eigenvalues and eigenvectors if not already done
        if ~isfield(solution, 'Uinv')
            Lmat0 = CreateSuperoperatorMatrix(@L0, input, solution);
            [U, D] = eig(full(Lmat0));
            Uinv = eye(M^2) / U;
        
            solution.U = U;
            solution.Uinv = Uinv;
            solution.D = sparse(D);
        end
        
        mat = LBTSS_Smat(input, solution);
        return
    end

    mat = zeros(M^2);
    for k = 1:M^2
        % create sampling matrices
        tmp = zeros(M^2, 1);
        tmp(k) = 1.0;
        tmp = reshape(tmp, [M M]);
        
        % enter the resulting matrix into mat in vector form
        mat(:, k) = reshape( L(input, tmp, solution), [M^2 1] );
    end

end