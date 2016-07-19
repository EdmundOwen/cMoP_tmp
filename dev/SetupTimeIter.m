%%% Setup the time iteration parameters and append them to the input 
%%% struct.

function input = SetupTimeIter(input, varargin)

    %% define default input values
    p = inputParser;
    
    p = AddParserOption(p, input, 'dt', 0.001, @isnumeric);
    p = AddParserOption(p, input, 'Nt', 1000, @isnumeric);
    p = AddParserOption(p, input, 'method', 'heun', true);
    p = AddParserOption(p, input, 'NTest', 1, @isnumeric);
    p = AddParserOption(p, input, 'probelist', {}, true);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.dt = p.Results.dt;
    input.Nt = p.Results.Nt;
    input.method = p.Results.method;
    input.NTest = p.Results.NTest;
    input.probelist = p.Results.probelist;
    
    %% add optional operators to probe
    
    %% if memory is needed, the probes must be made every time step
    if input.memoryNeeded
        input.probelist{end+1} = @SaveMemory;
        input.NTest = 1;
    end
    
    %% copy the time interval to the subinputs for each partition
    for k = 1:input.noPartitions
        input.subinput{k}.dt = input.dt;
    end
end
