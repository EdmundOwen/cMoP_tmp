%%% Test suite for the time iteration with the coupling to the environment
%%% appoximated to second order in the coupling, ie. including Born terms

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) TimeIterBTTest < matlab.unittest.TestCase
    
    properties
        rho
        input
        absTol = 1e-7;
        iterTol = 1e-3;     % tolerance for iterations
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2}};
        environ_varargin = {{}};
        timeiter_varargin = {{'operators', {@L0, @LMF, @LBT} }};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        
        function MethodSetup(tc, setup_varargin, environ_varargin, timeiter_varargin)
            tc.input.exists = true;
            tc.input = SetupSystem(tc.input, setup_varargin);
            tc.input = SetupEnvironment(tc.input, environ_varargin);
            tc.input = SetupTimeIter(tc.input, timeiter_varargin);
            tc.rho = InitializeRho(tc.input);
        end
        
    end
    
    methods (Test)
       
        function TestUnitaryEvolution(tc)
            %%% tests the unitary evolution using the model presented in
            %%% the original CMoP paper 
            %%% http://dx.doi.org/10.1103/PhysRevB.89.245108
            
            % setup the system as a unitary, undriven evolution of a spin
            % chain
            onsitedim = 2;
            Omega = 0.0;
            gamma = 0.0;
            U = 0.0;
            
            % create a set of inputs
            clustersize = 2;
            J = 0.5;
            Delta = 0.0;
            
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
            tc.rho = kron([0.0, 0.0; 0.0, 1.0], ...
                          [1.0, 0.0; 0.0, 0.0]);
            
            % do the iteration
            tc.input.dt = 0.01;
            tc.input.Nt = 600;
            result = TimeIter(tc.input, tc.rho);
            
            % test against the analytic solution by Flesch et al.
            % http://dx.doi.org/10.1103/PhysRevA.78.033608
            % Eq 15, f_{i,j}, i=j, i even
            exactocc = @(t) 0.5 + 0.5 * besselj(0, 4.0 * J * t);
            
            % calculate the occupation of the first site
            a_loc = annihilation(onsitedim);
            n1 = kron(a_loc' * a_loc, eye(onsitedim^(clustersize-1)));
            for i = 1:tc.input.Nt
                occ(i, 1) = abs(trace(n1 * result.hist{i}.rho));
                occ(i, 2) = exactocc(i * tc.input.dt);
                tc.assertEqual(occ(i, 1), occ(i, 2), 'AbsTol', tc.iterTol);
            end
            
        end
        
    end
    
end