close all; clear all;

addpath('../dev')
addpath('../utilities')

tic

%% physics
hbar = 1.0;
omega = 0.0;
Gamma_0 = 1.0;
k0 = (3.0 * pi * Gamma_0)^(1.0 / 3.0);
d_eg = 1.0;

%% scan variables
Delta_min = -1.0; nDelta = 41; deltaDelta = 0.05;
Delta_array = ((1:nDelta) - 1)*deltaDelta + Delta_min;

dr = 2.0
E = 0.1

%% Define spherical Bessel and Neumann functions
sphericalbesselj = @(nu, x) sqrt(pi ./ (2 .* x)) * besselj(nu + 0.5, x);
sphericalbessely = @(nu, x) sqrt(pi ./ (2 .* x)) * bessely(nu + 0.5, x);

%% parameters
% system
clustersize = 1;
onsitedim = 2;
L = { @L0, @LMF };
interaction_range = 1000;

% simulation
method = 'heun';
Nt = 3000;
dt = 0.02;
probelist = { @SaveDensity, @SaveMagnetisation, @SaveRho };
filebase = 'OpticalLattice_r%04.2f_E%03.1f.mat';

for Deltacount = 1:nDelta
    clear result
    clear input
    Delta = (Deltacount - 1) * deltaDelta + Delta_min;
    
    %% setup
    setup_varargin = {'clustersize', clustersize, 'onsitedim', onsitedim, 'L', L};
    environ_varargin = {};
    timeiter_varargin = { 'method', method, 'Nt', Nt, 'dt', dt, 'probelist', probelist };
    steady_varargin = {};

    input.exists = true;
    input = SetupSystem(input, setup_varargin);
    input = SetupEnvironment(input, environ_varargin);
    input = SetupTimeIter(input, timeiter_varargin);
    rho = InitializeRho(input);

    %% create time dependent Hamiltonian
    a = annihilation(onsitedim);
    % create local Hamiltonians
    H0_loc = 0.5 * hbar * Delta * (a * a' - a' * a);
    H0_mat = TensorAddOperator(H0_loc, onsitedim, clustersize);
    
    % create unitary couplings within the partition
    for i = 1:clustersize
        for j = 1:clustersize
            if i ~= j
                k0_r = k0 * abs(i - j) * dr;
                H0_mat = H0_mat + 0.5 * Gamma_0 * (sphericalbessely(0, k0_r) - 0.5 * sphericalbessely(2, k0_r)) ...
                                     * InsertOperator(a', i, onsitedim, clustersize) * InsertOperator(a, j, onsitedim, clustersize);
            end
        end
    end
    
    % create time-dependent pumping terms
    Evec = zeros(1, clustersize);
    for i = 1:clustersize
        Evec(i) = 0.5 * E;
    end
    Epumpup = zeros(onsitedim^clustersize);
    Epumpdown = zeros(onsitedim^clustersize);
    for i = 1:clustersize
        Epumpup = Epumpup + d_eg * Evec(i) * InsertOperator(a', i, onsitedim, clustersize);
        Epumpdown = Epumpdown + (d_eg * Evec(i))' * InsertOperator(a, i, onsitedim, clustersize);
    end
    input.H0 = @(t) H0_mat - exp(-1.0i * omega * t) * Epumpup - exp(1.0i * omega * t) * Epumpdown;

    input.subinput{1}.H0 = input.H0;

    %% and create Lindblad operators and weights
    A_Lindblad = cell(1, clustersize);
    Lindblad_weights = zeros(clustersize);
    for i = 1:clustersize
        for j = 1:clustersize
            A_Lindblad{i} = InsertOperator(a, i, onsitedim, clustersize);
            
            if i == j
                Lindblad_weights(i, j) = Gamma_0;
            else
                k0_r = k0 * abs(i - j) * dr;
                Lindblad_weights(i, j) = Gamma_0 * (sphericalbesselj(0, k0_r) - 0.5 * sphericalbesselj(2, k0_r));
            end
        end
    end
    % and allocate them
    input.subinput{1}.A_Lindblad = A_Lindblad;
    input.subinput{1}.Lindblad_weights = Lindblad_weights;

    %% create interactions
    int_type = {};
    keySet = {'interactionStrength', 'A', 'B', 'correlations', 'coordination', 'interactionType'};
    for i = 1:clustersize
        for j = 1:clustersize
            interaction_prefactor_unitary = 0.0;
            interaction_prefactor_dissipative = 0.0;
            
            % cycle through the interaction range checking whether this A/B
            % site pair is correct and adding to the interaction strength
            % to the associated prefactor
            for n = -interaction_range:interaction_range
                k0_r = k0 * abs(n) * dr;
            
                % continue if interaction is already modelled within the cluster
                if (i + n) > 0 && (i + n) < clustersize + 1
                    continue
                end
                % also continue if this isn't the interaction site we're
                % trying to find
                if j ~= mod(i + n, clustersize) + 1
                    continue
                end
            
                % add this interaction to the interaction prefactor
                interaction_prefactor_unitary = interaction_prefactor_unitary ...
                                                + 0.5 * Gamma_0 * (sphericalbessely(0, k0_r) - 0.5 * sphericalbessely(2, k0_r));
                interaction_prefactor_dissipative = interaction_prefactor_dissipative ...
                                                + Gamma_0 * (sphericalbesselj(0, k0_r) - 0.5 * sphericalbesselj(2, k0_r));
            end
            
            % create this interaction between A and B sites for all clusters
            B_SiteLabel = mod(i + j - 1, clustersize) + 1;
            correlations = ones(4*clustersize^2, 1);
            int_type{end+1} = containers.Map(keySet, {interaction_prefactor_dissipative, ...
                                                struct('Index', 1, 'SiteOperator', a, 'SiteLabel', i), ...
                                                struct('Index', 1, 'SiteOperator', a, 'SiteLabel', B_SiteLabel), ...
                                                correlations, 1, 'dissipative_ARB'});
            int_type{end+1} = containers.Map(keySet, {interaction_prefactor_dissipative, ...
                                                struct('Index', 1, 'SiteOperator', a, 'SiteLabel', i), ...
                                                struct('Index', 1, 'SiteOperator', a, 'SiteLabel', B_SiteLabel), ...
                                                correlations, 1, 'dissipative_BRA'});
            int_type{end+1} = containers.Map(keySet, {interaction_prefactor_unitary, ...
                                                struct('Index', 1, 'SiteOperator', a, 'SiteLabel', i), ...
                                                struct('Index', 1, 'SiteOperator', a', 'SiteLabel', B_SiteLabel), ...
                                                correlations, 1, 'unitary'});
            int_type{end+1} = containers.Map(keySet, {interaction_prefactor_unitary, ...
                                                struct('Index', 1, 'SiteOperator', a', 'SiteLabel', i), ...
                                                struct('Index', 1, 'SiteOperator', a, 'SiteLabel', B_SiteLabel), ...
                                                correlations, 1, 'unitary'});
        end
    end
    input.subinput{1}.interactions = CreateInteractions(input, int_type);

    %% run simulation
    % using steady state calculation
    solution = {};
    solution.rho = rho;
    input.H0 = input.H0(0);
    input.subinput{1}.H0 = input.subinput{1}.H0(0);
    input = SetupSteadyState(input, {});
    result = CalculateSteadyState(input, rho, solution);

    %% calculate cross-section
    sigma_tot = 0.0;
    % check that the iterator has converged, if not, set sigma to -2
    if ~isa(result, 'cell')
        if result == -1
            sigma_tot = -2.0 * clustersize;
        else
            throw exception
        end
    end
    for i = 1:clustersize
        sigma_tot = sigma_tot + imag(Evec(i) * d_eg ...
                            * trace(result{1} * InsertOperator(a, i, onsitedim, clustersize)));
    end
    sigma(Deltacount) = sigma_tot / clustersize;   
end

%% fit the cross-section to a Lorentzian to get half-decay rates and lineshifts for each set of lattice spacing
outputdata = lorentzfitMODIFIED(Delta_array, sigma, [], [], '3');
decayrate = 2.0*sqrt(outputdata{2}(3));
lineshift = outputdata{2}(2);

%% save data
filename = sprintf(filebase, dr, E);
save(filename, 'sigma', 'decayrate', 'lineshift');

toc
