%%% Setup the system.  This method defines the system, it's evolution
%%% operators and coordination... Basically all the possible inputs apart
%%% from the density matrix

function input = SetupSystem(input, varargin)

    %% define default input values
    p = inputParser;
    defaultClusterSize = 1;
    defaultOnsiteDim = 2;
    defaultDelta = 1.0;
    defaultJ = 1.0;
    defaultU = 0.0;
    defaultOmega = 0.0;
    defaultgamma = 0.0;
    defaultdim = 1;
    
    addOptional(p, 'clustersize', defaultClusterSize, @isnumeric);
    addOptional(p, 'onsitedim', defaultOnsiteDim, @isnumeric);
    addOptional(p, 'Delta', defaultDelta, @isnumeric);
    addOptional(p, 'J', defaultJ, @isnumeric);
    addOptional(p, 'U', defaultU, @isnumeric);
    addOptional(p, 'Omega', defaultOmega, @isnumeric);
    addOptional(p, 'gamma', defaultgamma, @isnumeric);
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
    
    %% initialise flag, the system is not a memory function so this is always false
    input.isMemory = false;
    
end