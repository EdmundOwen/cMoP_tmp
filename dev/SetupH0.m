%%% Setup a Hamiltonian.  This is a simple Hamiltonian which assumes
%%% an onsite energy term, a Hubbard U and a nearest neighbour coupling
%%% term

function H0 = SetupH0(input)

    %% get the useful values from the parser
    onsitedim = input.onsitedim;
    clustersize = input.clustersize;

    %% create a local on-site annihilation operator
    a_loc = annihilation(onsitedim);

    %% create a local on-site Hamiltonian
    Hloc = input.Delta * (a_loc')*a_loc ...
            + 0.5 * input.U * (a_loc')*(a_loc')*a_loc*a_loc ...
            + 0.5 * (input.Omega * (a_loc') + input.Omega' * a_loc);
    
    %% tensor product to give cluster onsite Hamiltonian
    H0 = kron(Hloc, speye(onsitedim^(clustersize-1)));
	for i = 2:clustersize
        H0 = H0 + kron(speye(onsitedim^(i - 1)), kron(Hloc, speye(onsitedim^(clustersize - i))));
    end
    
    %% add the couplings, which depends on the dimensionality of the system
    switch (input.dim)
        case 1
            Hintloc = -1.0 * input.J * (kron(a_loc', a_loc) + kron(a_loc, a_loc'));
            for i = 1:clustersize - 1
                H0 = H0 + kron(speye(onsitedim^(i - 1)), ...
                        kron(Hintloc, speye(onsitedim^(clustersize - i - 1))));
            end
            
        otherwise
            msgID = 'SetupH0:InvalidDimension';
            msg = sprintf('Error - Hamiltonians with dimension %i are not supported', input.dim);
            invalidinput_exception = MException(msgID, msg);
            throw(invalidinput_exception);
            
    end

end