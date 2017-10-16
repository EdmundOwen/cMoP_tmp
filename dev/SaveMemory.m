%%% A function to save the memory functions.  Only called if needed

function solution = SaveMemory( solution, input, timestep )
        
    for i = 1:input.noPartitions
        % get the density matrix
        %rho = solution.rho{i};
        subinput = input.subinput{i};
        partitionIndex = subinput.partitionIndex;

        for k = 1:numel(subinput.interactions)
            % get the interactions
            interaction_type = GetFromInput(subinput.interactions{k}, 'interactionType', 'unitary');
            Ak = subinput.interactions{k}.A; 
            Bk = subinput.interactions{k}.B;
            
            % get the density matrix on the partition which A and B are operating on
            rhoA = solution.rho{Ak.Index};
            rhoB = solution.rho{Bk.Index};

            % calculate the fluctuations
            dAk = Ak.Operator - trace(Ak.Operator * rhoA) * speye(size(Ak.Operator));
            dBk = Bk.Operator - trace(Bk.Operator * rhoB) * speye(size(Bk.Operator));

            % save the memory functions from this time step
            switch (interaction_type)
                case 'unitary'
                    solution.hist{partitionIndex}{timestep}.c{k} = dAk * rhoA;
                    solution.hist{partitionIndex}{timestep}.dkern{k} = dBk * rhoB;
                    solution.hist{partitionIndex}{timestep}.r{k} = rhoA * dAk;
                    solution.hist{partitionIndex}{timestep}.skern{k} = rhoB * dBk;
                case 'dissipative_ARB'
                    solution.hist{partitionIndex}{timestep}.c{k} = (Ak.Operator * rhoA - rhoA * Ak.Operator);
                    solution.hist{partitionIndex}{timestep}.dkern{k} = rhoB * dBk';
                    solution.hist{partitionIndex}{timestep}.r{k} = dAk * rhoA;
                    solution.hist{partitionIndex}{timestep}.skern{k} = (Bk.Operator' * rhoB - rhoB * Bk.Operator');
                case 'dissipative_BRA'
                    solution.hist{partitionIndex}{timestep}.c{k} = rhoA * dAk';
                    solution.hist{partitionIndex}{timestep}.dkern{k} = (Bk.Operator * rhoB - rhoB * Bk.Operator);
                    solution.hist{partitionIndex}{timestep}.r{k} = (Ak.Operator' * rhoA - rhoA * Ak.Operator');
                    solution.hist{partitionIndex}{timestep}.skern{k} = dBk * rhoB;
                otherwise
                    throw exception
            end                
        end
    end
    
end

