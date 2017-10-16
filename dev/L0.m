%%% Operation of the system evolution operator on the density matrix

function result = L0( input, mat, solution )

    %% if the input matrix is a struct, throw an error
    if isa(mat, 'struct')
        exception = MException('L0:InputInvalid', ...
            'the matrix input must not be a struct');
        throw(exception)
    end

    %% extract the Hamiltonian and Lindblad operators from the input
    H0 = input.H0;
    Lindblad_weights = input.Lindblad_weights;
    A_Lindblad = input.A_Lindblad;

    %% add unitary part of evolution to the superoperator
    if isa(full(H0), 'double')
        result = -1.0i * (H0 * mat - mat * H0);
    elseif isa(H0, 'function_handle')
        % assume that H0 is a time dependent function
        t = solution.time;
        result = -1.0i * (H0(t) * mat - mat * H0(t));
    else
        throw exception;
    end
    
    %% add dissipator
    [Ln_index, Lm_index, gammanm] = find(Lindblad_weights);
    for i = 1:numel(gammanm)
        Ln = A_Lindblad{Ln_index(i)};
        Lm = A_Lindblad{Lm_index(i)};
        result = result + 0.5 * gammanm(i) * (2.0 * Ln * mat * Lm' - mat * Lm' * Ln - Lm' * Ln * mat);
    end
    
end