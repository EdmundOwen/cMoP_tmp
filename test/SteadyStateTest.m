%%% Test suite for the steady state calculations with the coupling to the environment
%%% appoximated to second order in the coupling, ie. including Born terms

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) SteadyStateTest < matlab.unittest.TestCase
    
    properties
        input
        rho
        solution
        absTol = 1e-7;
        iterTol = 1e-4;
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 3, 'J', 0.5, 'Delta', 1.5, 'U', 10.0, 'Omega', 0.4, 'gamma', 0.3, 'operators', { @L0, @LMF } }, ...
                          {'gamma', 1.0, 'Omega', 1.0, 'Delta', 0.0, 'clustersize', 1, 'onsitedim', 2, 'operators', { @L0 } }, ...
                          {'Delta', 0.6, 'J', 1.0, 'gamma', 1.0, 'Omega', 1.5, 'clustersize', 2, 'onsitedim', 2, 'operators', { @L0, @LMF, @LBTSS }}};
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
            for k = 1:tc.input.noPartitions
            %%% create a free evolution superoperator matrix and test it
            M = tc.input.subinput{k}.M;
            
            for i = 1:numel(tc.input.subinput{k}.L)
                mat = CreateSuperoperatorMatrix( tc.input.subinput{k}.L{i}, tc.input.subinput{k}, tc.solution );
            
                % check that the matrices have the same size
                tc.assertEqual(size(mat), [M^2 M^2]);
                % and check that they do the right thing to rho
                tc.assertEqual(mat * reshape(tc.rho{k}, [M^2 1]), ...
                                reshape(tc.input.L{i}(tc.input.subinput{k}, tc.rho{k}, tc.solution), [M^2 1]), 'AbsTol', tc.absTol);
            end
            end
        end
        
        function TestSuperoperatorDiagonalization(tc)
            %%% tests the diagonalization of a superoperator matrix by
            %%% inputing the corresponding eigenvectors back into the
            %%% superoperators
            
            for k = 1:tc.input.noPartitions
            % setup
            M = tc.input.subinput{k}.M;
            
            for i = 1:numel(tc.input.L)
                % create a superoperator matrix and diagonalize
                Lmat = CreateSuperoperatorMatrix( tc.input.subinput{k}.L{i}, tc.input.subinput{k}, tc.solution );
                [U, D] = eig(full(Lmat));
                
                for j = 1:size(U, 1)
                    % input the jth eigenvector into the operator L
                    eigmat = reshape(U(:, j), [M M]);
                    result = tc.input.subinput{k}.L{i}(tc.input.subinput{k}, eigmat, tc.solution);
                    result = reshape(result, [M^2 1]);
                    
                    % this should be equal to the eigenvector multiplied by
                    % the eigenvalue, ie. L * U = D * U
                    tc.assertEqual(result, D(j, j) * U(:, j), 'AbsTol', tc.absTol);
                end
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
            
            % check to see that the inputs are valid
            if tc.input.onsitedim ~= 2 || tc.input.gamma == 0.0 || tc.input.Delta ~= 0.0
                return
            end
            
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
            for k = 1:tc.input.noPartitions
                tc.assertEqual(full(result{k}), expected, 'AbsTol', tc.absTol);
            end
            
        end
        
        function TestDefaults(tc)
            %%% a function to test the default settings to get a steady 
            %%% density matrix, no matter what they are.  this uses the 
            %%% time iteration schemes where the relaxation time taken from
            %%% the damping parameter
            
            % if using the steady state Born term superoperator, don't
            % bother
            if isequal(tc.input.L, { @L0, @LMF, @LBTSS })
                return
            end
            
            % check whether there is a steady state and don't perform this
            % test if there isn't
            if ~CheckForSteadyState(CreateSuperoperatorMatrix(@L0, tc.input, tc.solution));
                return
            end
            
            % setup the default iterator
            tc.input = SetupTimeIter(tc.input, {'method', 'heun'});
            tc.input.Nt = round(20.0 / (tc.input.gamma * tc.input.dt));
            
            % iterate the wavefunction
            result_t = TimeIter(tc.input, tc.rho);
            
            % calculate the steady state
            result_ss = CalculateSteadyState(tc.input, tc.rho, tc.solution);
            
            % compare the two results
            for k = 1:tc.input.noPartitions
                tc.assertEqual(full(result_ss{k}), result_t.rho{k}, 'AbsTol', tc.input.SSError);
            end
            
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
            
            % check to see that the inputs are valid ... note that this should correspond to
            % a unique solution part of the parameter space (see Fig. 1 of)
            % reference
            if tc.input.onsitedim < 3 || tc.input.J == 0.0 ...
                    || tc.input.Delta == 0.0 || tc.input.Omega == 0.0 ...
                    || tc.input.gamma == 0.0 || tc.input.U == 0.0
                return
            end
            
            % calculate the steady state solution directly using cMoP
            tc.solution.rho = CalculateSteadyState(tc.input, tc.rho, tc.solution);
            
            for k = 1:tc.input.noPartitions
                % get values needed to evaluate test
                onsitedim = tc.input.subinput{k}.onsitedim;
                clustersize = tc.input.subinput{k}.clustersize;
                gamma = tc.input.subinput{k}.gamma;
                U = tc.input.subinput{k}.U;
                Delta = tc.input.subinput{k}.Delta;
                J = tc.input.subinput{k}.J;
                Omega = tc.input.subinput{k}.Omega;
                
            % check that the values for the mean-field bosonic coherence s
            % and occupations are consistent with the theoretical values
            % for each site
            a_loc = annihilation(onsitedim);
            for i = 1:clustersize
                a_site = InsertOperator(a_loc, i, onsitedim, clustersize);

                % calculate the mean-field bosonic coherence
                a_mf = trace(a_site * tc.solution.rho{k});
                % and input it into the theoretical equation... it should
                % be zero
                c = -2.0 * (-1.0 * Delta + 0.5i * gamma) / U;
                % effective coupling to one site is multiplied by the
                % system coordination number
                Jeff = J * tc.input.coordination;  
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
                tc.assertEqual(trace(a_site' * a_site * tc.solution.rho{k}), a_occ, 'AbsTol', 10.0 * tc.input.SSError);
            end
            end
            
        end
        
        function TestSteadyStateBornTerm(tc)
            % tests the steady state solutions of the driven dissipative 
            % spin chain using the initial cMoP paper as a guide by testing
            % against results created by clusterSS
            
            % check that the input superoperators contains the Born terms
            % and that it's a spin chain
            if ~isequal(tc.input.L, { @L0, @LMF, @LBTSS }) || tc.input.onsitedim ~= 2
                return
            end
            
            % make sure that there is a steady state solution
            tc.assertTrue(CheckForSteadyState(CreateSuperoperatorMatrix(@L0, tc.input, tc.solution)));
            
            % load the results from Michael's code: clusterSS
            clusterSS = load('test_data/clusterSS_test.mat');
            
            % cycle through J
            for i = 1:clusterSS.NJ
                
                % assign a new J and reset the system
                varargin = {'clustersize', clusterSS.clustersize, ...
                            'onsitedim', clusterSS.onsitedim, ...
                            'J', (i-1) * clusterSS.dJ, ...
                            'Delta', clusterSS.Delta, ...
                            'U', clusterSS.U, ...
                            'Omega', clusterSS.Omega, ...
                            'gamma', clusterSS.gamma};
                tc.input = SetupSystem(tc.input, varargin);
                tc.input = SetupEnvironment(tc.input, {});
                
                % calculate the new steady state
                result = CalculateSteadyState(tc.input, tc.rho, tc.solution);
            
                
            for k = 1:tc.input.noPartitions
                % calculate the expectation of the annihilation operator 
                % for the first site
                a = annihilation(tc.input.subinput{k}.onsitedim);
                for j = 2:tc.input.subinput{k}.clustersize
                    a = kron(a, eye(tc.input.subinput{k}.onsitedim));
                end
                a_result = trace(a * result{k});
                
                % calculate the occupation of the first site
                n = a' * a;
                n_result = trace(n * result{k});
            
                % assert that these calculated values must be equal to
                % those obtained with clusterSS
                tc.assertEqual(a_result, clusterSS.a(i), 'AbsTol', tc.iterTol);
                tc.assertEqual(n_result, clusterSS.n(i), 'AbsTol', tc.iterTol);
            end
            
            end
            
        end
        
    end
    
end

