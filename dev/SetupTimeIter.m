%%% Setup the time iteration parameters and append them to the input 
%%% struct.

function input = SetupTimeIter(input, varargin)

    %% define default input values
    p = inputParser;
    defaultdt = GetFromInput(input, 'dt', 0.001);
    defaultNt = GetFromInput(input, 'Nt', 1000);
    defaultMethod = GetFromInput(input, 'method', 'crank-nicolson');
    
    addOptional(p, 'dt', defaultdt, @isnumeric);
    addOptional(p, 'Nt', defaultNt, @isnumeric);
    addOptional(p, 'method', defaultMethod);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.dt = p.Results.dt;
    input.Nt = p.Results.Nt;
    input.method = p.Results.method;
    
end
