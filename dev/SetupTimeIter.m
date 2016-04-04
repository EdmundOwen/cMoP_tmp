%%% Setup the time iteration parameters and append them to the input 
%%% struct.

function input = SetupTimeIter(input, varargin)

    %% define default input values
    p = inputParser;
    defaultdt = GetFromInput(input, 'dt', 0.001);
    defaultNt = GetFromInput(input, 'Nt', 1000);
    defaultMethod = GetFromInput(input, 'method', 'crank-nicolson');
    
    addOptional(p, 'dt', defaultdt, @isnumeric);
    addOptional(p, 'Nt', defaultNt, @isinteger);
    addOptional(p, 'method', defaultMethod);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.dt = p.Results.dt;
    input.Nt = p.Results.Nt;
    input.method = p.Results.method;

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
    
end
