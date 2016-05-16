%%% a function to add the interactions of the system with its environment
%%% into the input struct from a template defined by the interaction_list

function interactions = CreateInteractions( input, interaction_list )

    %% retrieve useful data from input
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    interactions = cell(numel(interaction_list), 1);

    %% create_interaction_operators for H_int_list
    for i = 1:numel(interaction_list)
        % extract interaction data
        A_site_label = interaction_list{i}{2};
        A_site_operator = interaction_list{i}{3};
        B_site_label = interaction_list{i}{4};
        B_site_operator = interaction_list{i}{5};    
    
        % create system interaction operator
        interaction_list{i}{7} = InsertOperator(A_site_operator, A_site_label, onsitedim, clustersize);
        % create environment interaction operator 
        % (assumes that the system is translationally invariant)   
        interaction_list{i}{8} = InsertOperator(B_site_operator, B_site_label, onsitedim, clustersize);
        
        % input into interactions struct
        interactions{i} = {interaction_list{i}{1}, interaction_list{i}{7}, ...
                                    interaction_list{i}{8}, interaction_list{i}{6}};
    end

end

