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
            for k = 1:tc.input.noPartitions
                tc.assertEqual(ndims(tc.rho{k}), 2);
            end
        end
        function TestTrace(tc)
            for k = 1:tc.input.noPartitions
                tc.assertEqual(trace(tc.rho{k}), 1.0, 'AbsTol', tc.absTol);
            end
        end
        function TestHermitian(tc)
            for k = 1:tc.input.noPartitions
                hermDiff = (tc.rho{k} - tc.rho{k}');
                tc.assertLessThan(max(abs(hermDiff(:))), tc.absTol);
            end
        end
        
        % test the system parameter setup
        function TestH0Dim(tc)
            for k = 1:tc.input.noPartitions
                % test to make sure that the system Hamiltonian has the same
                % dimension as the density matrix
                tc.assertEqual(size(tc.rho{k}), size(tc.input.subinput{k}.H0));
            end
        end
        function TestH0Hermitian(tc)
            for k = 1:tc.input.noPartitions
                % test to make sure that the system Hamiltonian is hermitian
                result = tc.input.subinput{k}.H0' - tc.input.subinput{k}.H0;
                tc.assertLessThan(max(abs(result(:))), tc.absTol);
            end
        end
        function TestLindbladNum(tc)
            for k = 1:tc.input.noPartitions
                % there should not be more than M^2-1 Lindblad operators
                % otherwise one of the operators is a linear combination of the
                % others and time is wasted
                tc.assertLessThanOrEqual(numel(tc.input.subinput{k}.A_Lindblad), tc.input.subinput{k}.M^2-1);
            end
        end
        function TestLindbladDim(tc)
            for k = 1:tc.input.noPartitions
                % test to make sure that all of the Lindblad operators have the
                % same dimension as the density matrix
                for i = 1:numel(tc.input.A_Lindblad)
                    tc.assertEqual(size(tc.rho{k}), size(tc.input.subinput{k}.A_Lindblad{i}));
                end
            end
        end
        function TestLindbladTrace(tc)
            for k = 1:tc.input.noPartitions
                % Lindblad operators should be traceless
                for i = 1:numel(tc.input.subinput{k}.A_Lindblad)
                    tc.assertEqual(trace(tc.input.subinput{k}.A_Lindblad{i}), 0.0, 'AbsTol', tc.absTol);
                end
            end
        end
        function TestLindbladWeights(tc)
            for k = 1:tc.input.noPartitions
                % the diagonalised Lindblad operator weights (damping 
                % factors) should be greater than or equal to zero
                tc.assertGreaterThanOrEqual(eig(tc.input.subinput{k}.Lindblad_weights), 0.0);
            end
        end
        
        % test parameter loading from file
        function TestLoad(tc)
            tc.input = LoadData('./test_data/test_load.mat', tc.input);
            tc.assertTrue(tc.input.flagLoadedData);
        end
        
    end
end