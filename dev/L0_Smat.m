%%% a function to create the superoperator matrix for the free evolution
%%% operator

function result = L0_Smat( input, solution )

    %% extract the Hamiltonian and Lindblad operators from the input
    M = input.M;
    H0 = input.H0;
    Lindblad_weights = input.Lindblad_weights;
    A_Lindblad = input.A_Lindblad;

    %% create prototypes for the superoperator matrix product functions
    % matprodleft and matprodright
    lintmp = reshape(1:M^2, [M M]);
    protoleft = kron(speye(M), lintmp);
    protoright = kron(lintmp', speye(M));
    
    %% add unitary part of evolution to the superoperator
    result = -1.0i * (matprod(H0, protoleft) - matprod(H0, protoright));
    
    %% add dissipator
    [Ln_index, Lm_index, gammanm] = find(Lindblad_weights);
    for i = 1:numel(gammanm)
        Ln = A_Lindblad{Ln_index(i)};
        Lm = A_Lindblad{Lm_index(i)};
        result = result + 0.5 * gammanm(i) * (2.0 * matprod(Ln, protoleft) * matprod(Lm', protoright) ...
                                                - matprod(Lm' * Ln, protoright) - matprod(Lm' * Ln, protoleft));
    end
    
end

