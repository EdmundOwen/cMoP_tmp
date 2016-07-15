%%% Test suite for the time iteration with mean-field coupling to the
%%% environment

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) TimeIterMFTest < matlab.unittest.TestCase
    
    properties
        rho
        input
        absTol = 1e-7;
        iterTol = 1e-3;     % tolerance for iterations
    end
    
    properties (MethodSetupParameter)
        setup_varargin = {{'clustersize', 1, 'onsitedim', 2, 'operators', { @L0, @LMF } },...
                          {'clustersize', 1, 'onsitedim', 3, 'J', 0.5, 'Delta', 1.5, 'U', 10.0, 'Omega', 0.4, 'gamma', 0.3, 'operators', { @L0, @LMF } }};
        environ_varargin = {{}};
        timeiter_varargin = {{'method', 'euler'},...
                             {'method', 'runge-kutta'},...
                             {'method', 'heun'}};
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
        
        function TestBoseHubbardMeanField(tc)
            % tests the mean-field, steady state solutions of the driven
            % dissipative Bose-Hubbard model using the results from Boite
            % et al. http://dx.doi.org/10.1103/PhysRevLett.110.233601
            
            % mean-field does not work for Crank-Nicolson so test for an
            % error in the test
            tc.assertNotEqual(tc.input.method, 'crank-nicolson');
            
            % check to see that the inputs are valid ... note that this should correspond to
            % a unique solution part of the parameter space (see Fig. 1 of)
            % reference
            if tc.input.onsitedim < 3 || tc.input.J == 0.0 ...
                    || tc.input.Delta == 0.0 || tc.input.Omega == 0.0 ...
                    || tc.input.gamma == 0.0 || tc.input.U == 0.0
                return
            end
            
            % do the iteration
            tc.input.Nt = round(20.0 / (tc.input.gamma * tc.input.dt));
            result = TimeIter(tc.input, tc.rho);
            
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
                a_mf = trace(a_site * result.rho{k});
                % and input it into the theoretical equation... it should
                % be zero
                c = -2.0 * (-1.0 * Delta + 0.5i * gamma) / U;
                % effective coupling to one site is multiplied by the
                % system coordination number
                Jeff = J * tc.input.subinput{k}.coordination;  
                % in the paper, the driving frequency term is not divided by 2
                F = 0.5 * Omega;
                amf_eq = @(a) ((F - Jeff * a) / (-1.0 * Delta + 0.5i * gamma) ...
                                * hypergeom([], [1+c, c'], 8.0 * abs((F - Jeff * a) / U)^2) ...
                                / hypergeom([], [c, c'], 8.0 * abs((F - Jeff * a) / U)^2)) - a;
                tc.assertEqual(amf_eq(a_mf), 0.0, 'AbsTol', tc.iterTol);
                
                % and the expected occupation (different from paper in the
                % hypergeometric function's third argument, but correct...)
                a_occ = abs(2.0 * (F - Jeff * a_mf) / U)^2 * ...
                            gammacomplex(c) * gammacomplex(c') / (gammacomplex(c+1) * gammacomplex(c'+1)) * ...
                            hypergeom([], [c+1, c'+1], 8.0 * abs((F - Jeff * a_mf) / U)^2) ...
                            / hypergeom([], [c, c'], 8.0 * abs((F - Jeff * a_mf) / U)^2);
                tc.assertEqual(trace(a_site' * a_site * result.rho{k}), a_occ, 'AbsTol', tc.iterTol);
            end
            end
            
        end
        
    end
    
end
