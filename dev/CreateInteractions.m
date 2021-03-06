%%% a function to add the interactions of the system with its environment
%%% into the input struct from a template defined by the interaction_list

function interactions = CreateInteractions( input, interaction_map )

    %% retrieve useful data from input
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    interactions = cell(numel(interaction_map), 1);

    %% create_interaction_operators for H_int_list
    for i = 1:numel(interaction_map)
        % extract interaction data
        A_site_label = interaction_map{i}('A').SiteLabel;
        A_site_operator = interaction_map{i}('A').SiteOperator;
        B_site_label = interaction_map{i}('B').SiteLabel;
        B_site_operator = interaction_map{i}('B').SiteOperator;
        
        % create system interaction operator
        A_operator = InsertOperator(A_site_operator, A_site_label, onsitedim, clustersize);
        % create environment interaction operator 
        % (assumes that the system is translationally invariant)   
        B_operator = InsertOperator(B_site_operator, B_site_label, onsitedim, clustersize);
        
        % input into interactions struct
        interactions{i}.interactionStrength = interaction_map{i}('interactionStrength');
        interactions{i}.A.Index = interaction_map{i}('A').Index;
        interactions{i}.A.Operator = A_operator;
        interactions{i}.B.Index = interaction_map{i}('B').Index;
        interactions{i}.B.Operator = B_operator;
        interactions{i}.Correlations = interaction_map{i}('correlations');
        
        % get the interaction type if it exists
        if isKey(interaction_map{i}, 'interactionType')
            interactions{i}.interactionType = interaction_map{i}('interactionType');
        else % default is unitary
            interactions{i}.interactionType = 'unitary';
        end
        
        % get the coordination number if it exists
        if isKey(interaction_map{i}, 'coordination')
            interactions{i}.coordination = interaction_map{i}('coordination');
        else
            interactions{i}.coordination = 1;
        end
    end

end

