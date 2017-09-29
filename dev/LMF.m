%%% Mean field superoperator

function result = LMF( input, mat, solution )

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
        interaction_strength = SE_interactions{i}.interactionStrength;
        A = SE_interactions{i}.A;
        B = SE_interactions{i}.B;
        interaction_type = GetFromInput(SE_interactions, 'interactionType', 'unitary');
        coordination = SE_interactions{i}.coordination;
    
        %% calculate the mean field term for this interaction
        switch (interaction_type)
            case 'unitary'
                result = result + -1.0i * interaction_strength * coordination * comm(B, A, mat, solution.rho);
            case 'dissipative'
                result = result + 0.5 * interaction_strength * coordination * (comm(B', A, mat, solution.rho) - comm(B, A', mat, solution.rho));
            otherwise
                throw exception
        end                
    
    end
    
end  % of function: LMF


%%% A sub-function to calculate the mean-field commutator, that is:
%%%
%%% Tr{B * rho} [A, mat]
function result = comm(B, A, mat, rho)

mean_field = trace(B.Operator * rho{B.Index});
result = mean_field * (A.Operator * mat - mat * A.Operator);

end  % of function: comm