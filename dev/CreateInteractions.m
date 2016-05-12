%%% a function to add the interactions of the system with its environment
%%% into the input struct from a template defined by the interaction_list

function input = CreateInteractions( input, interaction_list )

    %% retrieve useful data from input
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;

    %% create_interaction_operators for H_int_list
    for i = 1:numel(interaction_list)
        % extract interaction data
        A_site_label = interaction_list{i}{2};
        A_site_operator = interaction_list{i}{3};
        B_site_label = interaction_list{i}{4};
        B_site_operator = interaction_list{i}{5};    
    
        % create system interaction operator
        interaction_list{i}{7} = kron(speye(onsitedim^(A_site_label-1)), ...
                kron(A_site_operator, speye(onsitedim^(clustersize-A_site_label))));
        % create environment interaction operator 
        % (assumes that the system is translationally invariant)   
        interaction_list{i}{8} = kron(speye(onsitedim^(B_site_label - 1)), ...
                kron(B_site_operator, speye(onsitedim^(clustersize-B_site_label))));
        
        % input into interactions struct
        input.interactions{i} = {interaction_list{i}{1}, interaction_list{i}{7}, ...
                                    interaction_list{i}{8}, interaction_list{i}{6}};
    end

end

