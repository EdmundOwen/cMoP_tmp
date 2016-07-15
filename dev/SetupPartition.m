%%% Creates a sub input file which contains the information about the
%%% partition

function subinput = SetupPartition( input, partitionIndex )

    %% at the moment, all partitions are equal...
    subinput = input;
    subinput.partitionIndex = partitionIndex;

end

