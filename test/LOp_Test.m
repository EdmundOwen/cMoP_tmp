%%% Test suite for the system superoperator

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) LOp_Test < matlab.unittest.TestCase
    
    properties
        rho
        solution
        input
        absTol = 1e-7;
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2}};
        environ_varargin = {{}};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        
        function MethodSetup(tc, setup_varargin, environ_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.input = SetupEnvironment(tc.input, environ_varargin);
            tc.rho = InitializeRho(tc.input);
            tc.solution.rho = tc.rho;
        end
        
    end
    
    methods (Test)
        
        %% for L0
        % Check that the operator maintains the size of the density matrix
        function TestL0Size(tc)
            for k = 1:tc.input.noPartitions
                tc.assertEqual(size(L0(tc.input.subinput{k}, tc.rho{k}, tc.solution)), ...
                               size(tc.rho{k}));
            end
        end
        % Check that the operator result is Hermitian by operating on the set of 
        % unit vectors density matrices (ie. with zeros everywhere except 
        % a single 1 on the diagonal) and checking that the norm is 
        % preserved
        function TestL0Hermitian(tc)
            for i = 1:tc.input.M
                rho_unit = zeros(tc.input.M);
                rho_unit(i, i) = 1.0;
                tc.solution.rho = rho_unit;
                
                result = L0(tc.input, rho_unit, tc.solution);
                hermDiff = result' - result;
                tc.assertLessThan(max(abs(hermDiff(:))), tc.absTol);
            end
        end
        
        %% for LMF
        % Check that the operator maintains the size of the density matrix
        function TestLMFSize(tc)
            for k = 1:tc.input.noPartitions
                tc.assertEqual(size(LMF(tc.input.subinput{k}, tc.rho{k}, tc.solution)), ...
                               size(tc.rho{k}));
            end
        end
        % Check that the operator result is Hermitian by operating on the set of 
        % unit vectors density matrices (ie. with zeros everywhere except 
        % a single 1 on the diagonal) and checking that the norm is 
        % preserved
        function TestLMFHermitian(tc)
            for k = 1:tc.input.noPartitions
                for i = 1:tc.input.subinput{k}.M
                    rho_unit = zeros(tc.input.subinput{k}.M);
                    rho_unit(i, i) = 1.0;
                    tc.solution.rho{k} = rho_unit;
                
                    result = LMF(tc.input.subinput{k}, rho_unit, tc.solution);
                    hermDiff = result' - result;
                    tc.assertLessThan(max(abs(hermDiff(:))), tc.absTol);
                end
            end
        end
        
    end
    
end
        