%%% Setup the steady state search parameters and append them to the input 
%%% struct.

function input = SetupSteadyState(input, varargin)

    %% define default input values
    p = inputParser;
    
    p = AddParserOption(p, input, 'Niter', 100, @isnumeric);
    p = AddParserOption(p, input, 'SSError', 1e-5, @isnumeric);
    p = AddParserOption(p, input, 'searchCount', 5, @isnumeric);
    p = AddParserOption(p, input, 'searchInterval', [0.0 1.0], true);
    p = AddParserOption(p, input, 'optimalIncrementSearchMethod', 'binarychopExact', true);
    p = AddParserOption(p, input, 'steadyStateMixingRatio', 1.0, @isnumeric);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.Niter = round(p.Results.Niter);
    input.SSError = p.Results.SSError;
    input.searchCount = round(p.Results.searchCount);
    input.searchInterval = p.Results.searchInterval;
    input.optimalIncrementSearchMethod = p.Results.optimalIncrementSearchMethod;
    input.steadyStateMixingRatio = p.Results.steadyStateMixingRatio;
    
end
