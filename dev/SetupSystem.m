%%% Setup the system.  This method defines the system, it's evolution
%%% operators and coordination... Basically all the possible inputs apart
%%% from the density matrix

function input = SetupSystem(input, varargin)

    %% define default input values
    p = inputParser;
    defaultClusterSize = GetFromInput(input, 'clustersize', 1);
    defaultOnsiteDim = GetFromInput(input, 'onsitedim', 2);
    defaultCoordination = GetFromInput(input, 'coordination', 2);
    defaultDelta = GetFromInput(input, 'Delta', 1.0);
    defaultJ = GetFromInput(input, 'J', 1.0);
    defaultU = GetFromInput(input, 'U', 0.0);
    defaultOmega = GetFromInput(input, 'Omega', 0.0);
    defaultgamma = GetFromInput(input, 'gamma', 0.0);
    defaultOperators = GetFromInput(input, 'L', { @L0 });
    defaultdim = GetFromInput(input, 'dim', 1);
    
    addOptional(p, 'clustersize', defaultClusterSize, @isnumeric);
    addOptional(p, 'onsitedim', defaultOnsiteDim, @isnumeric);
    addOptional(p, 'coordination', defaultCoordination, @isnumeric);
    addOptional(p, 'Delta', defaultDelta, @isnumeric);
    addOptional(p, 'J', defaultJ, @isnumeric);
    addOptional(p, 'U', defaultU, @isnumeric);
    addOptional(p, 'Omega', defaultOmega, @isnumeric);
    addOptional(p, 'gamma', defaultgamma, @isnumeric);
    addOptional(p, 'operators', defaultOperators);
    addOptional(p, 'dim', defaultdim, @isnumeric);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover system parameters
    input.Delta = p.Results.Delta;
    input.J = p.Results.J;
    input.U = p.Results.U;
    input.Omega = p.Results.Omega;
    input.gamma = p.Results.gamma;
    input.dim = p.Results.dim;
    input.L = p.Results.operators;
   
    %% and the system coordination number
    input.coordination = p.Results.coordination;
    
    %% calculate the system Hilbert space dimension
    input.onsitedim = p.Results.onsitedim;
    input.clustersize = p.Results.clustersize;
    input.M = input.onsitedim ^ input.clustersize;
       
    %% create the system Hamiltonian
    input.H0 = SetupH0(input);
    
    %% create the system Lindblad operators
    input.Lindblad_weights = input.gamma * speye(input.clustersize);
    input.A_Lindblad = SetupLindblad(input);
    
    %% create load flag to show whether input contains loaded data
    input.flagLoadedData = false;
    
    %% create a list of the linear operators
    input.isLinear = @L0;
    
    %% check whether memory functions are needed
    input.memoryNeeded = false;
    input.termsNeedingMemory = { @LBT };
    for i = 1:numel(input.L)
        for j = 1:numel(input.termsNeedingMemory)
            if isequal(input.L{i}, input.termsNeedingMemory{j})
                input.memoryNeeded = true;
            end
        end
    end
    
    %% initialise flag, the system is not a memory function so this is always false
    input.isMemory = false;
    
end