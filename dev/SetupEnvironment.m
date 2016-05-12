%%% Setup the system-environment interactions.  This method defines the 
%%% coupling to the environment and the relevant system-environment
%%% bilinear operators

function input = SetupEnvironment( input, varargin )

    J = input.J;
    V = input.V;
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    
    % Create single site annihilation operator
    a_loc = annihilation(onsitedim);

%% Cluster interaction correlation array
% this details which interactions allowed to interact with each other.  for
% example, if a correlation occurs on the right hand side of a line
% segment, it cannot be transmitted back from the left hand side of the
% segment
% correlation = [1 1 0 0 0 0 0 0; ...
%                1 1 0 0 0 0 0 0; ...
%                0 0 1 1 0 0 0 0; ...
%                0 0 1 1 0 0 0 0; ...
%                0 0 0 0 1 1 0 0; ...
%                0 0 0 0 1 1 0 0; ...
%                0 0 0 0 0 0 1 1; ...
%                0 0 0 0 0 0 1 1];
           
correlation = [1 1 0 0; ...
               1 1 0 0; ...
               0 0 1 1; ...
               0 0 1 1];

%% Create a cell array containing the interactions with the environment
H_int_list = {};
H_int_list{1} = {-J, 1, a_loc, clustersize, a_loc', correlation(1, :)};
H_int_list{2} = {-J, 1, a_loc', clustersize, a_loc, correlation(2, :)};
H_int_list{3} = {-J, clustersize, a_loc', 1, a_loc, correlation(3, :)};
H_int_list{4} = {-J, clustersize, a_loc, 1, a_loc', correlation(4, :)};
% H_int_list{1} = {V, 1, a_loc, clustersize, a_loc', correlation(1, :)};
% H_int_list{2} = {V, 1, a_loc', clustersize, a_loc, correlation(2, :)};
% H_int_list{3} = {V, clustersize, a_loc', 1, a_loc, correlation(3, :)};
% H_int_list{4} = {V, clustersize, a_loc, 1, a_loc', correlation(4, :)};


    %% create_interaction_operators for H_int_list
    for i = 1:numel(H_int_list)
        % extract interaction data
        A_site_label = H_int_list{i}{2};
        A_site_operator = H_int_list{i}{3};
        B_site_label = H_int_list{i}{4};
        B_site_operator = H_int_list{i}{5};    
    
        % create system interaction operator
        H_int_list{i}{7} = kron(speye(onsitedim^(A_site_label-1)), ...
                kron(A_site_operator, speye(onsitedim^(clustersize-A_site_label))));
        % create environment interaction operator 
        % (assumes that the system is translationally invariant)   
        H_int_list{i}{8} = kron(speye(onsitedim^(B_site_label - 1)), ...
                kron(B_site_operator, speye(onsitedim^(clustersize-B_site_label))));
        
        % input into interactions struct
        input.interactions{i} = {H_int_list{i}{1}, H_int_list{i}{7}, H_int_list{i}{8}, H_int_list{i}{6}};
    end

end

