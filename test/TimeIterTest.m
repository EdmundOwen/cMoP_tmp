%%% Test suite for the time iteration (without coupling to environment)

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) TimeIterTest < matlab.unittest.TestCase
    
    properties
        rho
        input
        absTol = 1e-7;
        iterTol = 1e-3;     % tolerance for iterations
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2, 'L', { @L0 }}, ...
                          {'clustersize', 1, 'onsitedim', 2, 'L', { @L0 }, 'J', 0.0, 'Omega', 1.0, 'gamma', 0.0, 'Delta', 0.0}, ...
                          {'clustersize', 2, 'onsitedim', 2, 'L', { @L0 }, 'J', 1.0, 'Omega', 0.0, 'gamma', 0.0, 'Delta', 0.0}, ...
                          {'clustersize', 2, 'onsitedim', 2, 'L', { @L0 }, 'J', 0.0, 'Omega', 0.4, 'gamma', 1.5, 'Delta', 0.0}};
        timeiter_varargin = {{ 'method', 'euler' },...
                             ...%{'method', 'crank-nicolson'},...
                             {'method', 'runge-kutta'},...
                             {'method', 'heun'}};
    end

    methods (TestMethodSetup)
        
        function MethodSetup(tc, setup_varargin, timeiter_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.input = SetupTimeIter(tc.input, timeiter_varargin);
            tc.rho = InitializeRho(tc.input);
        end
        
    end
    
    methods (Test)
        
        function TestNoIter(tc)
            
            % test that function doesn't change the density matrix when Nt = 0
            tc.input.Nt = 0;
            results = TimeIter(tc.input, tc.rho);
            
            for k = 1:tc.input.noPartitions
               	diff = tc.rho{k} - results.rho{k};
                tc.assertEqual(max(abs(diff(:))), 0.0, 'AbsTol', tc.absTol);
            end
        end
        function TestIdentHamiltonian(tc)
        	% test that the density matrix doesn't change if the
            % Hamiltonian is set to the identity and all the Lindblad
            % weights are set to zero
            for k = 1:tc.input.noPartitions
                tc.input.subinput{k}.H0 = speye(size(tc.input.subinput{k}.H0));
                tc.input.subinput{k}.Lindblad_weights = zeros(size(tc.input.subinput{k}.Lindblad_weights));
            end
            
            results = TimeIter(tc.input, tc.rho);
            
            for k = 1:tc.input.noPartitions
                diff = tc.rho{k} - results.rho{k};
                tc.assertEqual(max(abs(diff(:))), 0.0, 'AbsTol', tc.absTol);
            end
        end
        
        function TestSimpleQubit(tc)
            % test a simple two-state, pure state driven evolution for a single
            % half-period of the expected oscillation
            
            % check to see that the inputs are valid
            if tc.input.onsitedim ~= 2 || tc.input.J ~= 0.0 ...
                    || tc.input.gamma ~= 0.0 || tc.input.Delta ~= 0.0 || tc.input.Omega == 0.0
                return
            end
            
            % start the qubit in the zero state
            expected = cell(1,tc.input.noPartitions);
            for k = 1:tc.input.noPartitions
                tc.rho{k} = 1;  expected{k} = 1;
                rho_zero = [1, 0; 0, 0];
                rho_one = [0, 0; 0, 1];
                for i = 1:tc.input.subinput{k}.clustersize
                    tc.rho{k} = kron(tc.rho{k}, rho_zero);
                    expected{k} = kron(expected{k}, rho_one);
                end
            end
            
            % iterate
            tc.input.Nt = round(pi / (tc.input.Omega * tc.input.dt));
            results = TimeIter(tc.input, tc.rho);
            
            for k = 1:tc.input.noPartitions
                % check that the qubit is in the one state
                tc.assertEqual(results.rho{k}, expected{k}, 'AbsTol', tc.iterTol);
            end
        end
        
        function TestCoupledQubit(tc)
            % test a simple two-state, evolution between two degenerate
            % coupled qubits for a single quarter-period of the expected 
            % oscillation so that the sites are entangled in the 
            % (01+10)/sqrt(2) state
            tc.input.onsitedim = 2;
            tc.input.clustersize = 2;
            
            % check to see that the inputs are valid
            if tc.input.onsitedim ~= 2 || tc.input.clustersize ~= 2 || tc.input.J == 0.0 ...
                    || tc.input.gamma ~= 0.0 || tc.input.Delta ~= 0.0 || tc.input.Omega ~= 0.0
                return
            end
            
            % start the qubit in the zero state
            rho_zero = [1, 0; 0, 0];
            rho_one = [0, 0; 0, 1];
            
            for k = 1:tc.input.noPartitions
                tc.rho{k} = kron(rho_one, rho_zero);
            end
            expected = [0,     0,    0,  0;
                        0,   0.5, 0.5i,  0;
                        0, -0.5i,  0.5,  0;
                        0,     0,    0,  0];
            
            % iterate
            tc.input.Nt = round(pi / (4.0 * tc.input.J * tc.input.dt));
            results = TimeIter(tc.input, tc.rho);
            
            for k = 1:tc.input.noPartitions
                % check that the qubits are in a maximally entangled state
                tc.assertEqual(results.rho{k}, expected, 'AbsTol', tc.iterTol);
            end
        end
        
        function TestDampedQubit(tc)
            % test a highly damped qubit to make sure the Lindblad
            % operators are working
            
            % check to see that the inputs are valid
            if tc.input.onsitedim ~= 2 || tc.input.J ~= 0.0 || tc.input.Delta ~= 0.0 || tc.input.gamma == 0.0
                return
            end
            
            % the expected final density matrix (from steady-state Lindblad
            % equation)
            expected = num2cell(ones(1, tc.input.noPartitions));
            rho00 = (tc.input.gamma^2 + tc.input.Omega^2) / (tc.input.gamma^2 + 2.0 * tc.input.Omega^2);
            rho11 = 1.0 - rho00;
            rho01 = -1.0i * tc.input.Omega / tc.input.gamma * (1.0 - 2.0 * rho00);
            rho10 = -1.0 * rho01;
            
            mixed = [rho00, rho01; rho10, rho11];
            for k = 1:tc.input.noPartitions
                for i = 1:tc.input.clustersize                
                    expected{k} = kron(expected{k}, mixed);
                end
            end
            
            % iterate
            tc.input.Nt = round(20.0 / (tc.input.gamma * tc.input.dt));
            results = TimeIter(tc.input, tc.rho);
            
            % check that the qubits are all close to a maximally mixed
            % state
            for k = 1:tc.input.noPartitions
                tc.assertEqual(results.rho{k}, expected{k}, 'AbsTol', tc.iterTol);
            end
        end
        
        % test to make sure that if the Hamiltonian changes, the time
        % evolution changes too
        function TestChangeH0(tc)
            for k = 1:tc.input.noPartitions
                % replace the Hamiltonian with the identity
                newH0 = speye(tc.input.subinput{k}.M);
                tc.input.subinput{k}.H0 = newH0;
                
                % and remove all Lindblad coefficients
                newLindbladWeights = zeros(size(tc.input.subinput{k}.Lindblad_weights));
                tc.input.subinput{k}.Lindblad_weights = newLindbladWeights;
            end
            
            % evolve using the identity
            results = TimeIter(tc.input, tc.rho);
            
            for k = 1:tc.input.noPartitions
                % check that the density matrix hasn't changed
                tc.assertEqual(abs(results.rho{k}), abs(tc.rho{k}), 'AbsTol', tc.iterTol);
            end
        end
    end
    
end

