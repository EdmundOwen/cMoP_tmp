%%% Setup the time iteration parameters and append them to the input 
%%% struct.

function input = SetupTimeIter(input, varargin)

    %% define default input values
    p = inputParser;
    defaultdt = GetFromInput(input, 'dt', 0.001);
    defaultNt = GetFromInput(input, 'Nt', 1000);
    defaultMethod = GetFromInput(input, 'method', 'crank-nicolson');
    defaultNTest = GetFromInput(input, 'NTest', 1);
    defaultProbes = GetFromInput(input, 'probes', {});
    
    addOptional(p, 'dt', defaultdt, @isnumeric);
    addOptional(p, 'Nt', defaultNt, @isnumeric);
    addOptional(p, 'method', defaultMethod);
    addOptional(p, 'NTest', defaultNTest, @isnumeric);
    addOptional(p, 'probes', defaultProbes);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.dt = p.Results.dt;
    input.Nt = p.Results.Nt;
    input.method = p.Results.method;
    input.NTest = p.Results.NTest;
    input.probelist = p.Results.probes;
    
    %% add optional operators to probe
    
    %% if memory is needed, the probes must be made every time step
    if input.memoryNeeded
        input.probelist{end+1} = @SaveMemory;
        input.NTest = 1;
    end
    
end
