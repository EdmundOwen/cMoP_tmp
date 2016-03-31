%%% a function to perform a time step on mat using the input superoperators
%
% INPUTS:
% mat           - the matrix to iterate
% solution      - the solution struct containing all memory functions, etc.
%
function [new_mat, solution] = PerformTimeStep( mat, solution, L, input, dt )

    method = input.method;
    
    %% switch between the various possible time-step methods
    switch method
        case 'euler'
            [new_mat, solution] = Euler( mat, solution, L, input, dt);
            
        case 'crank-nicolson'
            [new_mat, solution] = CrankNicolson( mat, solution, L, input, dt);
            
        case 'runge-kutta'
            [new_mat, solution] = rungekutta( mat, solution, L, input, dt);
            
        otherwise
            exception = MException('PerformTimeStep:InvalidTimeStepMethod', ...
                strcat('the method ', method, ' is not supported'));
            throw(exception)
    end

end

%%% a simple euler iteration scheme
function mat = rungakutta(mat, solution, L, input, dt)
    
    k0 = zeros(size(mat));
    k1 = zeros(size(mat));
    k2 = zeros(size(mat));
    k3 = zeros(size(mat));
        
    %% not implemented exception
    exception = MException('PerformTimeStep:NotImplementedException', ...
        strcat('Not implemented!'));
    throw(exception)
        
    % iterate the density matrix
    mat = mat + (dt / 6.0) * (k0 + 2.0 * k1 + 2.0 * k2 + k3);
    
end