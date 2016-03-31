%%% a function to perform a time step on mat using the input superoperators
%
% INPUTS:
% mat           - the matrix to iterate
% solution      - the solution struct containing all memory functions, etc.
%
function solution = PerformTimeStep( solution, L, input, dt )

    method = input.method;
    rho = solution.rho;
    
    %% switch between the various possible time-step methods
    switch method
        case 'euler'
            [new_rho, solution] = Euler( rho, solution, L, input, dt);
            
        case 'crank-nicolson'
            [new_rho, solution] = cranknicolson( rho, solution, L, input, dt);
            
        case 'runge-kutta'
            [new_rho, solution] = rungekutta( rho, solution, L, input, dt);
            
        otherwise
            exception = MException('PerformTimeStep:InvalidTimeStepMethod', ...
                strcat('the method ', method, ' is not supported'));
            throw(exception)
    end
    
    %% set the solution density matrix to its new value
    solution.rho = new_rho;

end

%%% a crank-nicolson scheme
function mat = cranknicolson(mat, solution, L, input, dt)

    % ensure that M is correct
    input.M = input.onsitedim^input.clustersize;

    % create superoperator matrices by putting matrices with one one into 
    % the superoperator method
    function Lmat = makeLmatrix( L )
        Lmat = zeros(input.M^2);
        for k = 1:input.M^2
            tmprho = zeros(input.M^2, 1);
            tmprho(k) = 1.0;
            tmprho = reshape(tmprho, [input.M input.M]);
            
            Lmat(:, k) = reshape( L(input, tmprho, solution), [input.M^2 1] );
        end
    end

    %% check if solution contains constant superoperator matrices and save
    if ~isfield(solution, 'Lmat')
        solution.Lmat = cell([numel(L) 1]);
    end
    for i = 1:numel(L)
        for j = 1:numel(input.termsNeedingMemory)
            if isempty(solution.Lmat{i}) && ~isequal(L{i}, input.termsNeedingMemory{j}) || ...
                                                        isequal(input.termsNeedingMemory{j}, L{i})
                % if the matrix doesn't exist and it doesn't need memory or if the 
                % term uses memory, save a version of it to solution
                solution.Lmat{i} = makeLmatrix( L{i} );
            end
        end
    end
    
    %% create the right- and left-hand side matrices of the crank nicolson method
    rhs = eye(input.M^2);
    lhs = eye(input.M^2);
    for i = 1:numel(L)
        rhs = rhs + 0.5 * dt * solution.Lmat{i};
        lhs = lhs - 0.5 * dt * solution.Lmat{i};
    end
    
    %% solve and reshape
    b = rhs * reshape(mat, [input.M^2 1]);
    a = lhs \ b;
    
    mat = reshape(a, [input.M input.M]);
    
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