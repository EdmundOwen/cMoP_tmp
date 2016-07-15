%%% Initializes the density matrix

function rho = InitializeRho(input)

    rho = cell(1,input.noPartitions);

    for k = 1:input.noPartitions
        rho{k} = zeros(input.subinput{k}.M);
        rho{k}(1, 1) = 1.0;
        rho{k}(2, 1) = -1.0i;
        rho{k}(1, 2) = 1.0i;
	end
    
end