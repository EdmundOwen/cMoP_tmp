%%% Test suite for the steady state calculations with the coupling to the environment
%%% appoximated to second order in the coupling, ie. including Born terms

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) SteadyStateTest < matlab.unittest.TestCase
    
    properties
        input
        rho
        absTol = 1e-7;
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2, 'operators', { @L0, @LMF } }};
        environ_varargin = {{}};
    end

    methods (TestMethodSetup)
        
        function MethodSetup(tc, setup_varargin, environ_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.input = SetupEnvironment(tc.input, environ_varargin);
            tc.rho = InitializeRho(tc.input);
        end
        
    end
    
    methods (Test)
        
        function TestSuperoperatorCreation(tc)
            % create a free evolution superoperator matrix
            solution = {};
            solution.rho = tc.rho;
            M = tc.input.M;
            
            for i = 1:numel(tc.input.L)
                mat = CreateSuperoperatorMatrix( tc.input.L{i}, tc.input, solution );
            
                % check that the matrices have the same size
                tc.assertEqual(size(mat), [M^2 M^2]);
                % and check that they do the right thing to rho
                tc.assertEqual(mat * reshape(tc.rho, [M^2 1]), ...
                                reshape(tc.input.L{i}(tc.input, tc.rho, solution), [M^2 1]));
            end
        end
        
    end
    
end

