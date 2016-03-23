%%% Setup a the Lindblad operators

function A_Lindblad = SetupLindblad(input)

    %% get the useful values from the parser
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;
    A_Lindblad = {};

	%% create single site annihilation operator
    a_loc = zeros(onsitedim);
    for i = 1:onsitedim-1
        a_loc(i, i+1) = sqrt(i);
    end
    
    %% tensor product to give cluster Lindblad operator
    for i = 1:clustersize
        A_Lindblad{i} = kron(speye(onsitedim^(i - 1)), kron(a_loc, speye(onsitedim^(clustersize - i))));
    end

end