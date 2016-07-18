%%% Test suite for the time iteration with the coupling to the environment
%%% appoximated to second order in the coupling, ie. including Born terms

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) TimeIterBTTest < matlab.unittest.TestCase
    
    properties
        rho
        input
        absTol = 1e-7;
        iterTol = 1e-2;     % tolerance for iterations
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2, 'L', { @L0, @LMF, @LBT } }, ...
                          {'clustersize', 2, 'onsitedim', 2, 'L', { @L0, @LMF, @LBT }, 'Omega', 0.0, 'gamma', 0.0, 'J', 1.0, 'Delta', 0.0 }};
        environ_varargin = {{}};
        timeiter_varargin = {{'method', 'euler', 'dt', 0.005 },...
                             {'method', 'runge-kutta', 'dt', 0.005},...
                             {'method', 'heun', 'dt', 0.005}};
    end

    methods (TestMethodSetup)
        
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
            
            % add SaveRho to the list of probes to evaluate the single-site
            % density
            tc.input.probelist{end+1} = @SaveRho;
            
            % check to see that the inputs are valid ..
            % the system is a unitary, undriven evolution of a spin
            % chain
            if tc.input.onsitedim ~= 2 || mod(tc.input.clustersize, 2) ~= 0 || tc.input.J == 0.0 ...
                    || tc.input.gamma ~= 0.0 || tc.input.Omega ~= 0.0
                return
            end
            
            % tensor product the two site density matrix to create the
            % system density matrix
            twosite_rho = kron([0.0, 0.0; 0.0, 1.0], ...
                          [1.0, 0.0; 0.0, 0.0]);
            for k = 1:tc.input.noPartitions
                tc.rho{k} = twosite_rho;
                for i = 1:(tc.input.subinput{k}.clustersize / 2) - 1
                    tc.rho{k} = kron(tc.rho{k}, twosite_rho);
                end
            end
            
            % do the iteration
            tc.input.Nt = 30;
            result = TimeIter(tc.input, tc.rho);
            
            % test against the analytic solution by Flesch et al.
            % http://dx.doi.org/10.1103/PhysRevA.78.033608
            % Eq 15, f_{i,j}, i=j, i even
            exactocc = @(t) 0.5 + 0.5 * besselj(0, 4.0 * tc.input.J * t);
            
            for k = 1:tc.input.noPartitions
                % calculate the occupation of the first site
                a_loc = annihilation(2);
                n1 = InsertOperator(a_loc' * a_loc, 1, 2, tc.input.subinput{k}.clustersize);
                for i = 1:tc.input.Nt
                    occ(i, 1) = abs(trace(n1 * result.hist{k}{i}.rho));
                    occ(i, 2) = exactocc(i * tc.input.dt);
                end
            
                tc.assertEqual(occ(i, 1), occ(i, 2), 'AbsTol', tc.iterTol);
            end
        end
        
    end
    
end