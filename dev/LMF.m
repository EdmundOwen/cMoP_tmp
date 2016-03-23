%%% Mean field superoperator

function result = LMF( input, rho )

    %% if the input rho is a struct, extract the density matrix
    if isa(rho, 'struct')
        rho = rho.rho;
    end
    
    result = zeros(size(rho));

    %% extract the interactions from the input
    SE_interactions = input.interactions;
    
    %% cycle through each of the terms in the interaction Hamiltonian
    for i = 1:numel(SE_interactions)
    
        %% extract interaction data
        interaction_strength = SE_interactions{i}{1};
        A = SE_interactions{i}{2};
        B = SE_interactions{i}{3};        
        
        %% calculate the mean field term for this interaction
        result = result + -1.0i * interaction_strength * comm(B, A, rho);
    
    end
    
end  % of function: LMF


%%% A sub-function to calculate the mean-field commutator, that is:
%%%
%%% Tr{B * rho} [A, rho]
function result = comm(B, A, rho)

mean_field = trace(B * rho);
result = mean_field * (A * rho - rho * A);

end  % of function: comm