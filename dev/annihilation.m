%%% a function to create an annihilation operator with the given onsite
%%% dimension

function result = annihilation( onsitedim )

    result = zeros(onsitedim);
    for j = 1:onsitedim-1
        result(j, j+1) = sqrt(j);
    end

end

