%%% Test suite for the time iteration (without coupling to environment)

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) TimeIterTest < matlab.unittest.TestCase
    
    properties
        rho
        input
        absTol = 1e-7;
        iterTol = 1e-3;     % tolerance for iterations
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2}};
        timeiter_varargin = {{'operators', { @L0 }, 'method', 'euler' },...
                             {'operators', { @L0 }, 'method', 'crank-nicolson'},...
                             {'operators', { @L0 }, 'method', 'runge-kutta'}};
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
           
           	diff = tc.rho - results.rho;
            tc.assertEqual(max(abs(diff(:))), 0.0, 'AbsTol', tc.absTol);
        end
        function TestIdentHamiltonian(tc)
        	% test that the density matrix doesn't change if the
            % Hamiltonian is set to the identity and all the Lindblad
            % weights are set to zero
            tc.input.H0 = speye(size(tc.input.H0));
            tc.input.Lindblad_weights = zeros(size(tc.input.Lindblad_weights));
            results = TimeIter(tc.input, tc.rho);
            
            diff = tc.rho - results.rho;
            tc.assertEqual(max(abs(diff(:))), 0.0, 'AbsTol', tc.absTol);
        end
        
        function TestSimpleQubit(tc)
            % test a simple two-state, pure state driven evolution for a single
            % half-period of the expected oscillation
            tc.input.onsitedim = 2;
            
            % set most of the inputs to zero
            tc.input.Delta = 0.0;
            tc.input.J = 0.0;
            tc.input.U = 0.0;
            tc.input.Lindblad_weights = zeros(size(tc.input.Lindblad_weights));
            
            % drive the qubit (if it's not already)
            if tc.input.Omega == 0.0
                tc.input.Omega = 1.0;
            end
            
            % setup the new Hamiltonian
            tc.input.H0 = SetupH0(tc.input);
            tc.input.Nt = round(pi / (tc.input.Omega * tc.input.dt));
            
            % start the qubit in the zero state
            tc.rho = 1;  expected = 1;
            rho_zero = [1, 0; 0, 0];
            rho_one = [0, 0; 0, 1];
            for i = 1:tc.input.clustersize
                tc.rho = kron(tc.rho, rho_zero);
                expected = kron(expected, rho_one);
            end
            
            % iterate
            results = TimeIter(tc.input, tc.rho);
            
            % check that the qubit is in the one state
            tc.assertEqual(results.rho, expected, 'AbsTol', tc.iterTol);
        end
        
        function TestCoupledQubit(tc)
            % test a simple two-state, evolution between two degenerate
            % coupled qubits for a single quarter-period of the expected 
            % oscillation so that the sites are entangled in the 
            % (01+10)/sqrt(2) state
            tc.input.onsitedim = 2;
            tc.input.clustersize = 2;
            
            % set most of the inputs to zero
            tc.input.Delta = 0.0;
            tc.input.U = 0.0;
            tc.input.Omega = 0.0;
            tc.input.Lindblad_weights = zeros(size(tc.input.Lindblad_weights));
            
            % couple the qubits
            if tc.input.J == 0.0
                tc.input.J = 1.0;
            end
            
            % setup the new Hamiltonian
            tc.input.H0 = SetupH0(tc.input);
            tc.input.Nt = round(pi / (4.0 * tc.input.J * tc.input.dt));
            
            % start the qubit in the zero state
            rho_zero = [1, 0; 0, 0];
            rho_one = [0, 0; 0, 1];
            
            tc.rho = kron(rho_one, rho_zero);
            expected = [0,     0,    0,  0;
                        0,   0.5, 0.5i,  0;
                        0, -0.5i,  0.5,  0;
                        0,     0,    0,  0];
            
            % iterate
            results = TimeIter(tc.input, tc.rho);
            
            % check that the qubits are in a maximally entangled state
            tc.assertEqual(results.rho, expected, 'AbsTol', tc.iterTol);
        end
        
        function TestDampedQubit(tc)
            % test a highly damped qubit to make sure the Lindblad
            % operators are working
            tc.input.onsitedim = 2;
            
            % set most of the inputs to zero
            tc.input.Delta = 0.0;
            tc.input.U = 0.0;
            tc.input.J = 0.0;
            
            % drive the qubits
            if tc.input.Omega == 0.0
                tc.input.Omega = 1.0;
            end
            
            % set the Lindblad weights
            if tc.input.gamma == 0.0
                tc.input.gamma = 1.0;
            end
            tc.input.Lindblad_weights = tc.input.gamma * eye(numel(tc.input.A_Lindblad));
            
            % setup the new Hamiltonian
            tc.input.H0 = SetupH0(tc.input);
            tc.input.Nt = round(10.0 / (tc.input.gamma * tc.input.dt));
            
            % the expected final density matrix (from steady-state Lindblad
            % equation)
            expected = 1;
            rho00 = (tc.input.gamma^2 + tc.input.Omega^2) / (tc.input.gamma^2 + 2.0 * tc.input.Omega^2);
            rho11 = 1.0 - rho00;
            rho01 = -1.0i * tc.input.Omega / tc.input.gamma * (1.0 - 2.0 * rho00);
            rho10 = -1.0 * rho01;
            
            mixed = [rho00, rho01; rho10, rho11];
            for i = 1:tc.input.clustersize
                expected = kron(expected, mixed);
            end
            
            % iterate
            results = TimeIter(tc.input, tc.rho);
            
            % check that the qubits are all close to a maximally mixed
            % state
            tc.assertEqual(results.rho, expected, 'AbsTol', tc.iterTol);
        end
        
    end
    
end

