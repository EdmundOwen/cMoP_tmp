%%% a quick method which checks whether there is a boolean in an input
%%% struct and returns that bool if it exists or returns false if it does
%%% not
function result = CheckBool( input, field, default )

    if isfield(input, field)
        if isa(input.(field), 'logical')
            result = input.(field);
            return
        end
    end
    
    result = default;

end

