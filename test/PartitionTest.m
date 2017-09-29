%%% Test suite for the partitioned calculations where a set of coupled cMoP
%%% partitions interact with each other (and potentially, an environment)

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) PartitionTest < matlab.unittest.TestCase
    
    properties
        input
        rho
        solution
        absTol = 1e-7;
        iterTol = 1e-2;
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'noPartitions', 2, 'L', { @L0, @LMF, @LBT } }, ...
                          {'noPartitions', 2, 'L', { @L0, @LMF, @LBT }, 'clustersize', 4}};
        steady_varargin = {{}};
    end

    methods (TestMethodSetup)
        
        function MethodSetup(tc, setup_varargin, steady_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.input = SetupSteadyState(tc.input, steady_varargin);
            
            % setup solution struct
            tc.solution = {};
        end
        
    end
    
    methods (Test)
        
        function TestCoupledQubitPartition(tc)
            % a test which couples two, undamped qubits in different cMoP
            % partitions
            
            % define and test some test specific inputs           
            if tc.input.onsitedim ~= 2 || tc.input.clustersize ~= 1 ...
                    || tc.input.Omega ~= 0.0 || tc.input.gamma ~= 0.0 ...
                    || tc.input.noPartitions ~= 2 
                return;
            end
            
            % perform a 2pi*3/4 rotation of the qubit
            dt = 0.01; Nt = round(3.0 * pi / (4.0 * tc.input.J * dt));
            timeiter_varargin = {'dt', dt, 'Nt', Nt, 'probelist', { @SaveDensity } };
            tc.input = SetupTimeIter(tc.input, timeiter_varargin);
            
            % set up the interactions between the partitions
            tc.input.subinput{1}.interactions = {};
            tc.input.subinput{2}.interactions = {};
            a = annihilation(tc.input.onsitedim);
            % for the first partition
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 1, 'Operator', a), 'B', struct('Index', 2, 'Operator', a'), 'Correlations', [1 1]);
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 1, 'Operator', a'), 'B', struct('Index', 2, 'Operator', a), 'Correlations', [1 1]);
            % for the second partition
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 2, 'Operator', a), 'B', struct('Index', 1, 'Operator', a'), 'Correlations', [1 1]);
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 2, 'Operator', a'), 'B', struct('Index', 1, 'Operator', a), 'Correlations', [1 1]);
            
            % set up the density matrix
            tc.rho = cell(1,2);
            tc.rho{1} = [0, 0; 0, 1]; tc.rho{2} = [1, 0; 0, 0];
            
            % iterate the density matrix such that the excitation moves
            % from one partition to the other
            results = TimeIter(tc.input, tc.rho);
            
            % test against the expected density
            for k = 1:tc.input.noPartitions
                tc.assertEqual(results.hist{k}{end}.n, 0.5, 'AbsTol', tc.iterTol);
            end
        end
        
        function TestBTUnitaryPartition(tc)
            % a test which calculates the unitary evolution after a quantum
            % quench of a set of two level systems using the results of
            % Flesch et al. but for two coupled, partitioned cMoP systems
            % http://dx.doi.org/10.1103/PhysRevB.89.245108
            
            % add SaveRho to the list of probes to evaluate the single-site
            % density
            dt = 0.005; Nt = 30;
            timeiter_varargin = {'dt', dt, 'Nt', Nt, 'probelist', { @SaveDensity } };
            tc.input = SetupTimeIter(tc.input, timeiter_varargin);
            
            % set up the couplings between the partitions
            tc.input.subinput{1}.interactions = {};
            tc.input.subinput{2}.interactions = {};
            a1 = kron(annihilation(tc.input.onsitedim), speye(tc.input.onsitedim^(tc.input.clustersize - 1)));
            am = kron(speye(tc.input.onsitedim^(tc.input.clustersize - 1)), annihilation(tc.input.onsitedim));
            % for the first partition
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 1, 'Operator', a1), 'B', struct('Index', 2, 'Operator', am'), 'Correlations', [1 1 1 1]);
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 1, 'Operator', a1'), 'B', struct('Index', 2, 'Operator', am), 'Correlations', [1 1 1 1]);
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 1, 'Operator', am), 'B', struct('Index', 2, 'Operator', a1'), 'Correlations', [1 1 1 1]);
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 1, 'Operator', am'), 'B', struct('Index', 2, 'Operator', a1), 'Correlations', [1 1 1 1]);
            % for the second partition
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 2, 'Operator', a1), 'B', struct('Index', 1, 'Operator', am'), 'Correlations', [1 1 1 1]);
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 2, 'Operator', a1'), 'B', struct('Index', 1, 'Operator', am), 'Correlations', [1 1 1 1]);
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 2, 'Operator', am), 'B', struct('Index', 1, 'Operator', a1'), 'Correlations', [1 1 1 1]);
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Index', 2, 'Operator', am'), 'B', struct('Index', 1, 'Operator', a1), 'Correlations', [1 1 1 1]);
            
            % check to see that the inputs are valid ..
            % the system is a unitary, undriven evolution of a spin
            % chain
            if tc.input.onsitedim ~= 2 || mod(tc.input.clustersize, 2) ~= 0 || tc.input.J == 0.0 ...
                    || tc.input.gamma ~= 0.0 || tc.input.Omega ~= 0.0
                return
            end
            
            % tensor product the two site density matrix to create the
            % system density matrix
            twosite_rho = kron([0.0, 0.0; 0.0, 1.0], ...
                          [1.0, 0.0; 0.0, 0.0]);
            for k = 1:tc.input.noPartitions
                tc.rho{k} = twosite_rho;
                for i = 1:(tc.input.subinput{k}.clustersize / 2) - 1
                    tc.rho{k} = kron(tc.rho{k}, twosite_rho);
                end
            end
            
            % do the iteration
            result = TimeIter(tc.input, tc.rho);
            
            % test against the analytic solution by Flesch et al.
            % http://dx.doi.org/10.1103/PhysRevA.78.033608
            % Eq 15, f_{i,j}, i=j, i even
            exactocc = @(t) 0.5 + 0.5 * besselj(0, 4.0 * tc.input.J * t);
            
            for k = 1:tc.input.noPartitions
                % calculate the occupation of the first site
                for i = 1:tc.input.Nt
                    occ(i, 1) = result.hist{k}{i}.n(1);
                    occ(i, 2) = exactocc(i * tc.input.dt);
                end
            
                tc.assertEqual(occ(i, 1), occ(i, 2), 'AbsTol', tc.iterTol);
            end
        end
        
    end
end

