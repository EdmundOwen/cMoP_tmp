function [n, a] = Test(  )

clear all
        environ_varargin = {};
        steady_varargin = {'Niter', round(1000), 'SSError', 1e-8};
        
    dJ = 0.1;
    NJ = 21;

    for i = 1:NJ
        
        setup_varargin = {'Delta', 0.6, 'J', (i-1)*dJ, 'gamma', 1.0, ...
            'Omega', 1.5, 'clustersize', 1, 'onsitedim', 2, 'operators', { @L0, @LMF, @LBTSS }};
            input.exists = true;
            
            input = SetupSystem(input, setup_varargin);
            input = SetupSystem(input, setup_varargin);
            input = SetupEnvironment(input, environ_varargin);
            input = SetupSteadyState(input, steady_varargin);
            rho = InitializeRho(input);
            
            % setup solution struct
            solution = {};
            solution.rho = rho;
            
            result = CalculateSteadyState(input, rho, solution);
            
            a_loc = annihilation(input.onsitedim);
%            a_loc = kron(a_loc, eye(2));
            n_loc = a_loc' * a_loc;
            
            n(i) = abs(trace(n_loc * result));
            a(i) = trace(a_loc * result);
            i
    end

    x = 0:dJ:dJ*(NJ-1);
    x = 2.0 * x;        % multiply by Z
    figure(1);
    plot(x, abs(n));
    figure(2);
    plot(x, real(a));
    figure(3);
    plot(x, imag(a));

end

