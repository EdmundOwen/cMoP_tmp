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

    %% actually create the interactions
    input.interactions = CreateInteractions(input, H_int_list);
    
    %% couple the partitions
    for j = 1:numel(input.interactions)
        input.interactions{j}.B.Index = 1;
    end

    %% and attach these interactions to each partition
    for k = 1:input.noPartitions
        input.subinput{k}.interactions = input.interactions;
    end
end

