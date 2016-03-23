%%% Test suite for the basic setup functions

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) SetupTest < matlab.unittest.TestCase
    
    properties
        rho
        input
        absTol = 1e-7;
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2}, ...
                          {'clustersize', 3, 'onsitedim', 3}};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        
        function MethodSetup(tc, setup_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.rho = InitializeRho(tc.input);
        end
        
    end
    
    methods (Test)
        
        % test the density matrix setup
        function TestDimension(tc)
            tc.assertEqual(ndims(tc.rho), 2);
        end
        function TestTrace(tc)
            tc.assertEqual(trace(tc.rho), 1.0, 'AbsTol', tc.absTol);
        end
        function TestHermitian(tc)
            hermDiff = (tc.rho - tc.rho');
            tc.assertLessThan(max(abs(hermDiff(:))), tc.absTol);
        end
        
        % test the system parameter setup
        function TestH0Dim(tc)
            % test to make sure that the system Hamiltonian has the same
            % dimension as the density matrix
            tc.assertEqual(size(tc.rho), size(tc.input.H0));
        end
        function TestH0Hermitian(tc)
            % test to make sure that the system Hamiltonian is hermitian
            result = tc.input.H0' - tc.input.H0;
            tc.assertLessThan(max(abs(result(:))), tc.absTol);
        end
        function TestLindbladNum(tc)
            % there should not be more than M^2-1 Lindblad operators
            % otherwise one of the operators is a linear combination of the
            % others and time is wasted
            tc.assertLessThanOrEqual(numel(tc.input.A_Lindblad), tc.input.M^2-1);
        end
        function TestLindbladDim(tc)
            % test to make sure that all of the Lindblad operators have the
            % same dimension as the density matrix
            for i = 1:numel(tc.input.A_Lindblad)
                tc.assertEqual(size(tc.rho), size(tc.input.A_Lindblad{i}));
            end
        end
        function TestLindbladTrace(tc)
            % Lindblad operators should be traceless
            for i = 1:numel(tc.input.A_Lindblad)
                tc.assertEqual(trace(tc.input.A_Lindblad{i}), 0.0, 'AbsTol', tc.absTol);
            end
        end
        function TestLindbladWeights(tc)
            % the Lindblad operator weights (damping factors) should be
            % greater than or equal to zero
            tc.assertGreaterThanOrEqual(tc.input.Lindblad_weights, 0.0);
        end
        
        % test parameter loading from file
        function TestLoad(tc)
            tc.input = LoadData('./test_data/test_load.mat', tc.input);
            tc.assertTrue(tc.input.flagLoadedData);
        end
        
    end
end