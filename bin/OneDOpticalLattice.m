close all; clear all;

addpath('../dev')

tic

%% physics
hbar = 1.0;
omega0 = 1.0;
omega = 3.0;
d_eg = 1.0;
E = 0.0;
Gamma_0 = 3.0;
k0 = 1.0;
dr = 100.0;

%% parameters
% system
clustersize = 1;
onsitedim = 2;
L = { @L0, @LMF, @LBT };
Delta = hbar * omega;
gamma = Gamma_0;
Omega = d_eg * E;

% simulation
method = 'heun';
Nt = 200;
dt = 0.02;
probelist = { @SaveDensity, @SaveMagnetisation, @SaveRho };

%% setup
setup_varargin = {'clustersize', clustersize, 'onsitedim', onsitedim, 'L', L, ...
                    'Delta', Delta, 'gamma', gamma, 'Omega', Omega};
environ_varargin = {};
timeiter_varargin = { 'method', method, 'Nt', Nt, 'dt', dt, 'probelist', probelist };
steady_varargin = {};

input.exists = true;
input = SetupSystem(input, setup_varargin);
input = SetupEnvironment(input, environ_varargin);
input = SetupTimeIter(input, timeiter_varargin);
rho = InitializeRho(input);

%% create time dependent Hamiltonian
energymat = -0.5 * hbar * omega0 * [1.0, 0.0; ...
                                    0.0, -1.0];
input.H0 = @(t) (energymat - [0.0, d_eg * E * exp(-1.0i * omega * t);...
                            (d_eg * E * exp(-1.0i * omega * t))', 0.0]);
input.subinput{1}.H0 = input.H0;

%% create interactions
int_type = {};
interaction_range = 1;
a = annihilation(onsitedim);
keySet = {'interactionStrength', 'A', 'B', 'correlations', 'coordination', 'interactionType'};
for i = 1:interaction_range
    k0_r = k0 * i * dr;
    int_type{end+1} = containers.Map(keySet, {Gamma_0 * (bessely(0, k0_r) - 0.5 * bessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        [1 1 1 1], 2, 'dissipative_ARB'});
	int_type{end+1} = containers.Map(keySet, {Gamma_0 * (bessely(0, k0_r) - 0.5 * bessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        [1 1 1 1], 2, 'dissipative_BRA'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (besselj(0, k0_r) - 0.5 * besselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', 1), ...
                                        [1 1 1 1], 2, 'unitary'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (besselj(0, k0_r) - 0.5 * besselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        [1 1 1 1], 2, 'unitary'});
end
input.subinput{1}.interactions = CreateInteractions(input, int_type);

%% run simulation
result = TimeIter(input, rho);

%% output results
time_array = (1:Nt)*dt;
for time = 1:Nt
    sigmaxfinal1(:, time) = result.hist{1}{time}.sigmax;
end
figure(1)
plot(time_array, sigmaxfinal1');

toc