%%% A function to save the density matrox.  Only called if needed

function result = SaveRho( result, input, i )
        
    % save the density matrix (for post-simulation analysis rather than
    % calculations)
    rho = result.rho;
    
    % save this to the history
    result.hist{i}.rho = rho;

end

