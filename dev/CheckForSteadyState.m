%%% Check whether Lmat0 is a relaxing map a la Rivas and Huelga in Open 
%%% Quantum Systems - An Introduction, (2011), that is: for a relaxing map,
%%% the zero eigenvalue of the free evolution must not be degenerate and 
%%% all other eigenvalues must have negative real components

function [steadystate, errortype] = CheckForSteadyState( Lmat )

    %% machine tolerance error for calculating the zeros of the eigenvalue 
    % decomposition
    machine_tol = 1e-10;

    % default is true
    steadystate = true;
    errortype = 'none';

    % calculate the eigenvalues of the input matrix and test them
    eigvals = eig(full(Lmat));
    if sum(abs(eigvals) < machine_tol) > 1
        steadystate = false;
        errortype = 'zero eigenvalue is degenerate';
    elseif sum(abs(eigvals) < machine_tol) == 0
        steadystate = false;
        errortype = 'no zero eigenvalue';
    elseif sum(real(eigvals) > machine_tol) > 0
        steadystate = false;
        errortype = sprintf('%i eigenvalues have a positive real part', sum(real(eigvals) > machine_tol));
    end
    
end

