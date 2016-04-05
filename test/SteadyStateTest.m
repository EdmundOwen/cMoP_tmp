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
        setup_varargin = {{'gamma', 1.0, 'Omega', 1.0, 'clustersize', 4, 'onsitedim', 2, 'operators', { @L0 } }, ...
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
        
        function TestDefaults(tc)
            %%% a function to test the default settings to get a steady 
            %%% density matrix, no matter what they are.  this uses the 
            %%% time iteration schemes where the relaxation time is assumed
            %%% to be twenty times the maximum energy difference
            
            % check whether there is a steady state and don't perform this
            % test if there isn't
            if ~CheckForSteadyState(CreateSuperoperatorMatrix(@L0, tc.input, tc.solution));
                return
            end
            
            % calculate the spectrum of the isolated system
            energies = eig(tc.input.H0);
            
            % setup the default iterator
            tc.input = SetupTimeIter(tc.input, {'method', 'heun'});
            tc.input.Nt = 20 * round(abs(energies(1) - energies(end))) / tc.input.dt;
            
            % iterate the wavefunction
            result_t = TimeIter(tc.input, tc.rho);
            
            % calculate the steady state
            result_ss = CalculateSteadyState(tc.input, tc.rho, tc.solution);
            
            % compare the two results
            tc.assertEqual(result_ss, result_t.rho, 'AbsTol', tc.input.SSError);
            
        end
        
        function TestBoseHubbardMeanField(tc)
            % tests the mean-field, steady state solutions of the driven
            % dissipative Bose-Hubbard model using the results from Boite
            % et al. http://dx.doi.org/10.1103/PhysRevLett.110.233601
            
            % only do this test if the test case is in the mean field
            % approximation
            if ~isequal(tc.input.L, { @L0, @LMF })
                return
            end
            
            % create a set of inputs... note that this should correspond to
            % a unique solution part of the parameter space (see Fig. 1 of)
            % reference
            clustersize = 1;
            onsitedim = 3;     % caution: this needs to be greater than 2 to pass
            J = 0.5;
            Delta = 1.5;
            U = 10.0;
            Omega = 0.4;
            gamma = 0.3;
            
            % reset the system
            varargin = {'clustersize', clustersize, ...
                        'onsitedim', onsitedim, ...
                        'J', J, ...
                        'Delta', Delta, ...
                        'U', U, ...
                        'Omega', Omega, ...
                        'gamma', gamma};
            tc.input = SetupSystem(tc.input, varargin);
            tc.input = SetupEnvironment(tc.input, {});
            tc.rho = InitializeRho(tc.input);
            
            % calculate the steady state solution directly using cMoP
            tc.solution.rho = CalculateSteadyState(tc.input, tc.rho, tc.solution);
            
            % check that the values for the mean-field bosonic coherence s
            % and occupations are consistent with the theoretical values
            % for each site
            a_loc = annihilation(onsitedim);
            for i = 1:clustersize
                a_site = kron(speye(onsitedim^(i - 1)), ...
                            kron(a_loc, speye(onsitedim^(clustersize-i))));

                % calculate the mean-field bosonic coherence
                a_mf = trace(a_site * tc.solution.rho);
                % and input it into the theoretical equation... it should
                % be zero
                c = -2.0 * (-1.0 * Delta + 0.5i * gamma) / U;
                % effective coupling to one site (factor of 0.5 as only the couplings into the system count for MF)
                Jeff = J * (0.5 * numel(tc.input.interactions));  
                % in the paper, the driving frequency term is not divided by 2
                F = 0.5 * Omega;
                amf_eq = @(a) ((F - Jeff * a) / (-1.0 * Delta + 0.5i * gamma) ...
                                * hypergeom([], [1+c, c'], 8.0 * abs((F - Jeff * a) / U)^2) ...
                                / hypergeom([], [c, c'], 8.0 * abs((F - Jeff * a) / U)^2)) - a;
                % make the tolerance an order of magnitude larger than the
                % steady state error
                tc.assertEqual(amf_eq(a_mf), 0.0, 'AbsTol', 10.0 * tc.input.SSError);
                
                % and the expected occupation (different from paper in the
                % hypergeometric function's third argument, but correct...)
                a_occ = abs(2.0 * (F - Jeff * a_mf) / U)^2 * ...
                            gammacomplex(c) * gammacomplex(c') / (gammacomplex(c+1) * gammacomplex(c'+1)) * ...
                            hypergeom([], [c+1, c'+1], 8.0 * abs((F - Jeff * a_mf) / U)^2) ...
                            / hypergeom([], [c, c'], 8.0 * abs((F - Jeff * a_mf) / U)^2);
                % and similarly with the tolerance for the occupation 
                % expectation value
                tc.assertEqual(trace(a_site' * a_site * tc.solution.rho), a_occ, 'AbsTol', 10.0 * tc.input.SSError);
            end
            
        end
        
    end
    
end

