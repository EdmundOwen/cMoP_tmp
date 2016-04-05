%%% Setup the steady state search parameters and append them to the input 
%%% struct.

function input = SetupSteadyState(input, varargin)

    %% define default input values
    p = inputParser;
    defaultNiter = GetFromInput(input, 'Niter', 100);
    defaultSSError = GetFromInput(input, 'SSError', 1e-5);
    
    addOptional(p, 'Niter', defaultNiter, @isnumeric);
    addOptional(p, 'SSError', defaultSSError, @isnumeric);
    
    %% parse the inputs
    tmp = varargin{:};
    parse(p, tmp{:});
    
    %% recover iteration parameters
    input.Niter = round(p.Results.Niter);
    input.SSError = p.Results.SSError;
    
end
