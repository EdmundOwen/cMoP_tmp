%%% adds a parser option

function p = AddParserOption( p, input, optionname, defaultvalue, checkfunction )

    %% if no check is allowed (or valid), create a trivial checkfunction
    if checkfunction == true
        checkfunction = @(value) true;
    end

    valuefrominput = GetFromInput(input, optionname, defaultvalue);
    addOptional(p, optionname, valuefrominput, checkfunction);

end