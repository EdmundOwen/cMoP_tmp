close all; clear all;

addpath('../dev')

tic

%% physics
hbar = 1.0;
omega = 0.0;
d_eg = 1.0;
E = 0.1;
Gamma_0 = 2.0;
k0 = 1.0;
dr = 1000.0;
Delta = -0.5;

%% parameters
% system
clustersize = 2;
onsitedim = 2;
L = { @L0, @LMF };

% simulation
method = 'heun';
Nt = 600;
dt = 0.02;
probelist = { @SaveDensity, @SaveMagnetisation, @SaveRho };

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
H0_loc = -0.5 * hbar * Delta * (a * a' - a' * a);
H0_mat = TensorAddOperator(H0_loc, onsitedim, clustersize);
% create unitary couplings within the partition
if clustersize == 2
    H0_mat = H0_mat + 0.5 * Gamma_0 * (besselj(0, k0 * dr) - 0.5 * besselj(2, k0 * dr)) * (kron(a, a') + kron(a', a));
elseif clustersize > 2
    throw exception;
end
% create time-dependent pumping terms
Evec = zeros(1, clustersize);
for j = 1:clustersize
    Evec(j) = (-1.0)^j * E;
end
Epumpup = zeros(onsitedim^clustersize);
Epumpdown = zeros(onsitedim^clustersize);
for i = 1:clustersize
    Epumpup = Epumpup + d_eg * Evec(i) * InsertOperator(a', i, onsitedim, clustersize);
    Epumpdown = Epumpdown + (d_eg * Evec(i))' * InsertOperator(a, i, onsitedim, clustersize);
end
input.H0 = @(t) H0_mat - exp(-1.0i * omega * t) * Epumpup - exp(1.0i * omega * t) * Epumpdown;

input.subinput{1}.H0 = input.H0;

%% and create Lindblad operators
A_Lindblad = cell(1, clustersize);
Lindblad_weights = zeros(clustersize);
for i = 1:clustersize-1
    A_Lindblad{i} = InsertOperator(a, i, onsitedim, clustersize);
    Lindblad_weights(i, i) = Gamma_0;
    Lindblad_weights(i+1, i) = Gamma_0 * (bessely(0, k0 * dr) - 0.5 * bessely(2, k0 * dr));
    Lindblad_weights(i, i+1) = Gamma_0 * (bessely(0, k0 * dr) - 0.5 * bessely(2, k0 * dr));
end
A_Lindblad{clustersize} = InsertOperator(a, clustersize, onsitedim, clustersize);
Lindblad_weights(clustersize, clustersize) = Gamma_0;
% and allocate them
input.subinput{1}.A_Lindblad = A_Lindblad;
input.subinput{1}.Lindblad_weights = Lindblad_weights;

%% create interactions
int_type = {};
interaction_range = 1;
a = annihilation(onsitedim);
keySet = {'interactionStrength', 'A', 'B', 'correlations', 'coordination', 'interactionType'};
for i = 1:interaction_range
    k0_r = k0 * i * dr;
    int_type{end+1} = containers.Map(keySet, {Gamma_0 * (bessely(0, k0_r) - 0.5 * bessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        [1 1 1 1 0 0 0 0], 1, 'dissipative_ARB'});
	int_type{end+1} = containers.Map(keySet, {Gamma_0 * (bessely(0, k0_r) - 0.5 * bessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        [1 1 1 1 0 0 0 0], 1, 'dissipative_BRA'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (besselj(0, k0_r) - 0.5 * besselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', clustersize), ...
                                        [1 1 1 1 0 0 0 0], 1, 'unitary'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (besselj(0, k0_r) - 0.5 * besselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        [1 1 1 1 0 0 0 0], 1, 'unitary'});
    int_type{end+1} = containers.Map(keySet, {Gamma_0 * (bessely(0, k0_r) - 0.5 * bessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        [0 0 0 0 1 1 1 1], 1, 'dissipative_ARB'});
	int_type{end+1} = containers.Map(keySet, {Gamma_0 * (bessely(0, k0_r) - 0.5 * bessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        [0 0 0 0 1 1 1 1], 1, 'dissipative_BRA'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (besselj(0, k0_r) - 0.5 * besselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', 1), ...
                                        [0 0 0 0 1 1 1 1], 1, 'unitary'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (besselj(0, k0_r) - 0.5 * besselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', clustersize), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        [0 0 0 0 1 1 1 1], 1, 'unitary'});
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