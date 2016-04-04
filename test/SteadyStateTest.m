%%% Test suite for the steady state calculations with the coupling to the environment
%%% appoximated to second order in the coupling, ie. including Born terms

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) SteadyStateTest < matlab.unittest.TestCase
    
    properties
        input
        rho
        solution
        absTol = 1e-7;
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 4, 'onsitedim', 2, 'operators', { @L0 } }, ...
                          {'clustersize', 1, 'onsitedim', 2, 'operators', { @L0, @LMF } }};
        environ_varargin = {{}};
        steady_varargin = {{}};
    end

    methods (TestMethodSetup)
        
        function MethodSetup(tc, setup_varargin, environ_varargin, steady_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.input = SetupEnvironment(tc.input, environ_varargin);
            tc.input = SetupSteadyState(tc.input, steady_varargin);
            tc.rho = InitializeRho(tc.input);
            
            % setup solution struct
            tc.solution = {};
            tc.solution.rho = tc.rho;
        end
        
    end
    
    methods (Test)
        
        function TestSuperoperatorCreation(tc)
            %%% create a free evolution superoperator matrix and test it
            M = tc.input.M;
            
            for i = 1:numel(tc.input.L)
                mat = CreateSuperoperatorMatrix( tc.input.L{i}, tc.input, tc.solution );
            
                % check that the matrices have the same size
                tc.assertEqual(size(mat), [M^2 M^2]);
                % and check that they do the right thing to rho
                tc.assertEqual(mat * reshape(tc.rho, [M^2 1]), ...
                                reshape(tc.input.L{i}(tc.input, tc.rho, tc.solution), [M^2 1]));
            end
        end
        
        function TestSuperoperatorDiagonalization(tc)
            %%% tests the diagonalization of a superoperator matrix by
            %%% inputing the corresponding eigenvectors back into the
            %%% superoperators
            
            % setup
            M = tc.input.M;
            
            for i = 1:numel(tc.input.L)
                % create a superoperator matrix and diagonalize
                Lmat = CreateSuperoperatorMatrix( tc.input.L{i}, tc.input, tc.solution );
                [U, D] = eig(Lmat);
                
                for j = 1:size(U, 1)
                    % input the jth eigenvector into the operator L
                    eigmat = reshape(U(:, j), [M M]);
                    result = tc.input.L{i}(tc.input, eigmat, tc.solution);
                    result = reshape(result, [M^2 1]);
                    
                    % this should be equal to the eigenvector multiplied by
                    % the eigenvalue, ie. L * U = D * U
                    tc.assertEqual(result, D(j, j) * U(:, j), 'AbsTol', tc.absTol);
                end
            end
            
        end
        
        function TestDampedQubit(tc)
            %%% tests the steady state density matrix for the damped qubit
            %%% as in TimeIterTest
            
            % test for setups containing only linear terms
            if ~isequal(tc.input.L, { @L0 })
                return
            end
            
            % it's a qubit...
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
            
            % calculate the steady state solution
            result = CalculateSteadyState(tc.input, tc.rho, tc.solution);
            
            % test to make sure that the steady state solution is correct
            tc.assertEqual(expected, result, 'AbsTol', tc.absTol);
            
        end
        
    end
    
end

