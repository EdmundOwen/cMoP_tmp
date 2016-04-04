%%% a function to create superoperator matrices by putting matrices with a 
%%% single one into the superoperator method

function mat = CreateSuperoperatorMatrix( L, input, solution )

    mat = zeros(input.M^2);
    for k = 1:input.M^2
        % create sampling matrices
        tmp = zeros(input.M^2, 1);
        tmp(k) = 1.0;
        tmp = reshape(tmp, [input.M input.M]);
        
        % enter the resulting matrix into mat in vector form
        mat(:, k) = reshape( L(input, tmp, solution), [input.M^2 1] );
	end

end

