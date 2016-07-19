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
        setup_varargin = {{'noPartitions', 2, 'L', { @L0, @LMF, @LBT } }};
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
                    || tc.input.noPartitions ~= 2
                return;
            end
            
            dt = 0.01; Nt = round(pi / (4.0 * tc.input.J * dt));
            timeiter_varargin = {'dt', dt, 'Nt', Nt, 'probelist', { @SaveDensity } };
            tc.input = SetupTimeIter(tc.input, timeiter_varargin);
            
            % set up the interactions between the partitions
            tc.input.subinput{1}.interactions = {};
            tc.input.subinput{2}.interactions = {};
            a = annihilation(tc.input.onsitedim);
            % for the first partition
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Operator', a), 'B', struct('Index', 2, 'Operator', a'), 'Correlations', [1 1]);
            tc.input.subinput{1}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Operator', a'), 'B', struct('Index', 2, 'Operator', a), 'Correlations', [1 1]);
            % for the second partition
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Operator', a), 'B', struct('Index', 1, 'Operator', a'), 'Correlations', [1 1]);
            tc.input.subinput{2}.interactions{end+1} ...
                = struct('interactionStrength', tc.input.J, 'A', struct('Operator', a'), 'B', struct('Index', 1, 'Operator', a), 'Correlations', [1 1]);
            
            % set up the density matrix
            tc.rho = cell(1,2);
            tc.rho{1} = [0, 0; 0, 1]; tc.rho{2} = [1, 0; 0, 0];
            
            % iterate the density matrix such that the excitation moves
            % from one partition to the other
            results = TimeIter(tc.input, tc.rho);
            
            %%tmp
            for i = 1:tc.input.Nt
                occ(i) = trace(results.hist{1}{i}.n);
            end
            
            % test against the expected density
            a = annihilation(tc.input.onsitedim);
            for k = 1:tc.input.noPartitions
                tc.assertEqual(results.hist{k}{end}.n, 0.5, 'AbsTol', tc.iterTol);
            end
        end
        
        function TestCoupledDampedQubitPartition(tc)
            % a test which finds the steady state for two, damped qubits 
            % in different cMoP partitions
        end
        
        function TestPartitionvsNoPartition(tc)
            % a test which calculates the density distribution from the
            % default parameters and compares the result when there is
            % partitioning against when there is none
        end
        
        function TestBTUnitaryPartition(tc)
            % a test which calculates the unitary evolution after a quantum
            % quench of a set of two level systems using the results of
            % Flesch et al. but for two coupled, partitioned cMoP systems
        end
    end
end

