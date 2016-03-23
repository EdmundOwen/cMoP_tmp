%%% A function which loads data from a file and puts it into the input 
%%% struct.  All data loaded should have been created using the SetupSystem
%%% method so we check this.  The input contains a flag for loaded data
%%% which is switched when this method is run

function input = LoadData( filename, input )

    %% load data and get filenames
    datain = load(filename);
    labels = fieldnames(datain);
    
    for i = 1:numel(labels)
        tmp_label = labels(i);
        tmp_label = tmp_label{1};
        
        %% check that the label is valid by comparing it to the input
        if isfield(input, tmp_label)
            input.(tmp_label) = datain.(tmp_label);
        else
            msgID = 'LoadData:InvalidInputLabel';
            msg = strcat('Error - ', filename, ' contains invalid input data');
            invalidinput_exception = MException(msgID, msg);
            throw(invalidinput_exception);
        end
    end
    
    input.flagLoadedData = true;

end

