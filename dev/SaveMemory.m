%%% A function to save the memory functions.  Only called if needed

function solution = SaveMemory( solution, input, timestep )
        
    for i = 1:input.noPartitions
        % get the density matrix
        %rho = solution.rho{i};
        subinput = input.subinput{i};
        partitionIndex = subinput.partitionIndex;

        for k = 1:numel(subinput.interactions)
            % get the interactions
            Ak = subinput.interactions{k}.A; 
            Bk = subinput.interactions{k}.B;
            
            % get the density matrix on the partition which A and B are operating on
            rhoA = solution.rho{Ak.Index};
            rhoB = solution.rho{Bk.Index};

            % calculate the fluctuations
            dAk = Ak.Operator - trace(Ak.Operator * rhoA) * speye(size(Ak.Operator));
            dBk = Bk.Operator - trace(Bk.Operator * rhoB) * speye(size(Bk.Operator));

            % save the memory functions from this time step
            solution.hist{partitionIndex}{timestep}.c{k} = dAk * rhoA;
            solution.hist{partitionIndex}{timestep}.dkern{k} = dBk * rhoB;
            solution.hist{partitionIndex}{timestep}.r{k} = rhoA * dAk;
            solution.hist{partitionIndex}{timestep}.skern{k} = rhoB * dBk;
        end
    end
    
end

