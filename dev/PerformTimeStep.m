%%% a function to perform a time step on mat using the input superoperators
%
% INPUTS:
% mat           - the matrix to iterate
% solution      - the solution struct containing all memory functions, etc.
%
function mat = PerformTimeStep(mat, solution, L, input )

    method = input.method;
    
    %% switch between the various possible time-step methods
    switch method
        case 'euler'
            mat = euler(mat, solution, L, input);
            
        case 'runge-kutta'
            mat = rungekutta(mat, solution, L, input);
            
        otherwise
            exception = MException('PerformTimeStep:InvalidTimeStepMethod', ...
                strcat('the method ', method, ' is not supported'));
            throw(exception)
            
    end

end

%%% a simple euler iteration scheme
function mat = euler(mat, solution, L, input)
    
    rhs = zeros(size(mat));
    % calculate rhs of Louiville equaion
    for j = 1:numel(L)
        rhs = rhs + L{j}( input, mat, solution );
    end
        
    % iterate the density matrix
    mat = mat + input.dt * rhs;
    
end

%%% a crank-nicolson scheme
function mat = cranknicolson(mat, solution, L, input)

    %% not implemented exception
    exception = MException('PerformTimeStep:NotImplementedException', ...
        strcat('Not implemented!'));
    throw(exception)
        
end

%%% a simple euler iteration scheme
function mat = rungakutta(mat, solution, L, input)
    
    k0 = zeros(size(mat));
    k1 = zeros(size(mat));
    k2 = zeros(size(mat));
    k3 = zeros(size(mat));
        
    %% not implemented exception
    exception = MException('PerformTimeStep:NotImplementedException', ...
        strcat('Not implemented!'));
    throw(exception)
        
    % iterate the density matrix
    mat = mat + (input.dt / 6.0) * (k0 + 2.0 * k1 + 2.0 * k2 + k3);
    
end