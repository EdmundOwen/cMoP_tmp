function rho_array = Test(  )

clear all

    addpath('../dev');
    
        environ_varargin = {};
        steady_varargin = {'Niter', round(1000), 'SSError', 1e-8};
        
    dV = 0.2;
    NV = 1;
    
    dU = 0.2;
    NU = 1;
    
    rho_array = {};

    for i = 1:NV
        for j = 1:NU
        
        setup_varargin = {'Delta', 0.0, 'J', 0.0, 'gamma', 1.0, 'U', (j-1)*dU, 'V', (i-1)*dV, ...
            'Omega', 0.8, 'clustersize', 2, 'onsitedim', 3, 'operators', { @L0, @LMF, @LBTSS }};
            input.exists = true;
            
            input = SetupSystem(input, setup_varargin);
            input = SetupSystem(input, setup_varargin);
            input = SetupEnvironment(input, environ_varargin);
            input = SetupSteadyState(input, steady_varargin);
            if i == 1 && j == 1
                rho = InitializeRho(input);
            else
                rho = solution.rho;
            end
            
            % setup solution struct
            solution = {};
            solution.rho = rho;
            
            result = CalculateSteadyState(input, rho, solution);
            
            rho_array{i, j} = result;
            
            strcat(sprintf('%i', i), ',', sprintf('%i', j));
        end
    end

%     x = 0:dV:dV*(NV-1);
%     x = 2.0 * x;        % multiply by Z
%     figure(1);
%     plot(x, abs(n));
%     figure(2);
%     plot(x, real(a));
%     figure(3);
%     plot(x, imag(a));

end

