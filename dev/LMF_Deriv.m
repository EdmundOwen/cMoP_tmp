%%% Jacobian of the mean field superoperator matrix

function result = LMF_Deriv( input, mat, solution )

    %% if the input matrix is a struct, throw an error
    if isa(mat, 'struct')
        exception = MException('L0:InputInvalid', ...
            'the matrix input must not be a struct');
        throw(exception)
    end
    
    if input.noPartitions ~= 1
        exception = MException('LMF_Deriv:NotImplemented', ...
            'derivative methods for coupled partitions not implemented');
        throw(exception)
    end
    
    M = input.M;
    result = zeros(M^2);
    
    %% extract the interactions from the input
    SE_interactions = input.interactions;
    
    %% cycle through each of the terms in the interaction Hamiltonian
    for i = 1:numel(SE_interactions)
    
        %% reset the result for this interaction
        iterresult = zeros(M^2);
        
        %% extract interaction data
        interactionStrength = SE_interactions{i}.interactionStrength;
        A = SE_interactions{i}.A;
        B = SE_interactions{i}.B;
        interaction_type = GetFromInput(SE_interactions, 'interactionType', 'unitary');
    
        %% cycle through projection matrices
        for k = 1:M^2
            % create sampling matrices
            Mk = zeros(M^2, 1);
            Mk(k) = 1.0;
            Mk = reshape(Mk, [M M]);
            
            tmp = -1i * interactionStrength * ...
                        (trace(B.Operator * Mk) * (A.Operator * mat - mat * A.Operator) ...
                        + trace(B.Operator * mat) * (A.Operator * Mk - Mk * A.Operator));
                    
            %% reshape and input into result
            iterresult(:, k) = reshape(tmp, [M^2 1]);
        end
        
        result = result + iterresult;
    
    end
    
end  % of function: LMF_Deriv