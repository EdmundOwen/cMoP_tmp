%%% inserts an operator into a lattice at the given position

function result = InsertOperator( operator, position, onsitedim, clustersize )

    result = kron(speye(onsitedim^(position - 1)), ...
              kron(operator, speye(onsitedim^(clustersize - position))));

end

