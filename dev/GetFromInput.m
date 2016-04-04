%%% a function to get values from the input with possible default value
function value = GetFromInput(input, name, default)

    if isfield(input, name)
        value = input.(name);
    else
        value = default;
    end
        
end