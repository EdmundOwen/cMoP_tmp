%%% Setup the system.  This method defines the system, it's evolution
%%% operators and coordination... Basically all the possible inputs apart
%%% from the density matrix

function input = SetupSystem(input, varargin)

    %% define default input values
    p = inputParser;
    
    %% add the parser options and defaults
    p = AddParserOption(p, input, 'clustersize', 1, @isnumeric);
    p = AddParserOption(p, input, 'onsitedim', 2, @isnumeric);
    p = AddParserOption(p, input, 'coordination', 2, @isnumeric);
    p = AddParserOption(p, input, 'Delta', 1.0, @isnumeric);
    p = AddParserOption(p, input, 'J', 1.0, @isnumeric);
    p = AddParserOption(p, input, 'U', 0.0, @isnumeric);
    p = AddParserOption(p, input, 'V', 0.0, @isnumeric);
    p = AddParserOption(p, input, 'Omega', 0.0, @isnumeric);
    p = AddParserOption(p, input, 'gamma', 0.0, @isnumeric);
    p = AddParserOption(p, input, 'L', { @L0 }, true);
    p = AddParserOption(p, input, 'dim', 1, @isnumeric);
    p = AddParserOption(p, input, 'noPartitions', 1, @isnumeric);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover system parameters
    input.Delta = p.Results.Delta;
    input.J = p.Results.J;
    input.U = p.Results.U;
    input.V = p.Results.V;
    input.Omega = p.Results.Omega;
    input.gamma = p.Results.gamma;
    input.dim = p.Results.dim;
    input.L = p.Results.L;
    input.noPartitions = p.Results.noPartitions;
   
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
            if strcmp(func2str(input.L{i}), func2str(input.termsNeedingMemory{j}))
                input.memoryNeeded = true;
            end
        end
    end
    
    %% initialise flag, the system is not a memory function so this is always false
    input.isMemory = false;
    
    %% set up the partitions using this data
    for i = 1:input.noPartitions
        input.subinput{i} = SetupPartition(input, i);
    end
end