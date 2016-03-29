%%% Setup the time iteration parameters and append them to the input 
%%% struct.

function input = SetupTimeIter(input, varargin)

    %% define default input values
    p = inputParser;
    defaultdt = 0.001;
    defaultNt = 1000;
    defaultOperators = { @L0 };
    defaultMethod = 'euler';
    
    addOptional(p, 'dt', defaultdt, @isnumeric);
    addOptional(p, 'Nt', defaultNt, @isinteger);
    addOptional(p, 'operators', defaultOperators);
    addOptional(p, 'method', defaultMethod);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.dt = p.Results.dt;
    input.Nt = p.Results.Nt;
    input.L = p.Results.operators;
    input.method = p.Results.method;

    %% check whether memory functions are needed
    input.memoryNeeded = false;
    termsNeedingMemory = { @LBT };
    for i = 1:numel(input.L)
        for j = 1:numel(termsNeedingMemory)
            if isequal(input.L{i}, termsNeedingMemory{j})
                input.memoryNeeded = true;
            end
        end
    end
    
end
