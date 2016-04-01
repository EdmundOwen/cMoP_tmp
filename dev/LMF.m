%%% Mean field superoperator

function result = LMF( input, mat, solution )

    %% if the input rho is a struct, extract the density matrix
    %% if the input matrix is a struct, throw an error
    if isa(mat, 'struct')
        exception = MException('L0:InputInvalid', ...
            'the matrix input must not be a struct');
        throw(exception)
    end
    
    result = zeros(size(mat));

    %% extract the interactions from the input
    SE_interactions = input.interactions;
    
    %% cycle through each of the terms in the interaction Hamiltonian
    for i = 1:numel(SE_interactions)
    
        %% extract interaction data
        interaction_strength = SE_interactions{i}{1};
        A = SE_interactions{i}{2};
        B = SE_interactions{i}{3};        
        
        %% calculate the mean field term for this interaction
        result = result + -1.0i * interaction_strength * comm(B, A, mat, solution.rho);
    
    end
    
end  % of function: LMF


%%% A sub-function to calculate the mean-field commutator, that is:
%%%
%%% Tr{B * rho} [A, mat]
function result = comm(B, A, mat, rho)

mean_field = trace(B * rho);
result = mean_field * (A * mat - mat * A);

end  % of function: comm