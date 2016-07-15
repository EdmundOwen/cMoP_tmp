%%% a function to perform a time step on mat using the input superoperators
%
% INPUTS:
% mat           - the matrix to iterate
% solution      - the solution struct containing all memory functions, etc.
%
function [new_mat, solution] = PerformTimeStep( mat, solution, L, input, dt, method )
    
    %% switch between the various possible time-step methods
    switch method
        case 'euler'
            [new_mat, solution] = Euler( mat, solution, L, input, dt);
            
        case 'crank-nicolson'
            [new_mat, solution] = CrankNicolson( mat, solution, L, input, dt);
            
        case 'runge-kutta'
            [new_mat, solution] = RungeKutta( mat, solution, L, input, dt);
            
        case 'heun'
            [new_mat, solution] = Heun( mat, solution, L, input, dt);
            
        otherwise
            exception = MException('PerformTimeStep:InvalidTimeStepMethod', ...
                strcat('the method ', method, ' is not supported'));
            throw(exception)
    end

end