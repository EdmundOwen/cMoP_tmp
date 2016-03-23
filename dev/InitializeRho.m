%%% Initializes the density matrix

function rho = InitializeRho(input)

    rho = zeros(input.M);
    rho(1, 1) = 1.0;
    rho(2, 1) = -1.0i;
    rho(1, 2) = 1.0i;
    
end