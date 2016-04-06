function result = matprod( mat, tmp )

    M = max(size(mat));
    [i, j, val] = find(tmp);
    
    result = sparse(i, j, mat(val), M^2, M^2);
    
end