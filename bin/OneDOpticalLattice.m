close all; clear all;

addpath('../dev')
addpath('../utilities')

tic

%% physics
hbar = 1.0;
omega = 0.0;
E = 0.1;
Gamma_0 = 1.0;
k0 = 2.0 * pi;
d_eg = 1.0;
% variables
Delta_min = -1.1; Delta_start = 1; nDelta = 21; deltaDelta = 0.1;
r_min = 0.0; r_start = 1; nr = 22; deltar = 0.2;

%% Define spherical Bessel and Neumann functions
sphericalbesselj = @(nu, x) sqrt(pi ./ (2 .* x)) * besselj(nu + 0.5, x);
sphericalbessely = @(nu, x) sqrt(pi ./ (2 .* x)) * bessely(nu + 0.5, x);

%% parameters
% system
clustersize = 1;
onsitedim = 2;
L = { @L0, @LMF };

% simulation
method = 'heun';
Nt = 3000;
dt = 0.02;
probelist = { @SaveDensity, @SaveMagnetisation, @SaveRho };

for rcount = r_start:nr
    dr = rcount * deltar + r_min
for Deltacount = Delta_start:nDelta
    clear result
    clear input
    Delta = Deltacount * deltaDelta + Delta_min;
    
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
if clustersize == 2
    H0_mat = H0_mat + 0.5 * Gamma_0 * (sphericalbessely(0, k0 * dr) ...
                                        - 0.5 * sphericalbessely(2, k0 * dr)) * (kron(a, a') + kron(a', a));
elseif clustersize > 2
    throw exception;
end
% create time-dependent pumping terms
Evec = zeros(1, clustersize);
for j = 1:clustersize
    Evec(j) = E;
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
    Lindblad_weights(i+1, i) = Gamma_0 * (sphericalbesselj(0, k0 * dr) - 0.5 * sphericalbesselj(2, k0 * dr));
    Lindblad_weights(i, i+1) = Gamma_0 * (sphericalbesselj(0, k0 * dr) - 0.5 * sphericalbesselj(2, k0 * dr));
end
A_Lindblad{clustersize} = InsertOperator(a, clustersize, onsitedim, clustersize);
Lindblad_weights(clustersize, clustersize) = Gamma_0;
% and allocate them
input.subinput{1}.A_Lindblad = A_Lindblad;
input.subinput{1}.Lindblad_weights = Lindblad_weights;

%% create interactions
int_type = {};
interaction_range = 3;
a = annihilation(onsitedim);
keySet = {'interactionStrength', 'A', 'B', 'correlations', 'coordination', 'interactionType'};
% create interaction dictionary
for i = 1:interaction_range
    k0_r = k0 * i * dr;
    % create correlation vectors
    correlations_Left = horzcat(zeros(1, 4*(i-1)), ones(1, 4), zeros(1, 4), zeros(1, 4*(interaction_range-i)));
    correlations_Right = horzcat(zeros(1, 4*(i-1)), zeros(1, 4), ones(1, 4), zeros(1, 4*(interaction_range-i)));
    % create individual interactions
    int_type{end+1} = containers.Map(keySet, {Gamma_0 * (sphericalbesselj(0, k0_r) - 0.5 * sphericalbesselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        correlations_Left, 2, 'dissipative_ARB'});
	int_type{end+1} = containers.Map(keySet, {Gamma_0 * (sphericalbesselj(0, k0_r) - 0.5 * sphericalbesselj(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        correlations_Left, 2, 'dissipative_BRA'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (sphericalbessely(0, k0_r) - 0.5 * sphericalbessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', clustersize), ...
                                        correlations_Left, 2, 'unitary'});
    int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (sphericalbessely(0, k0_r) - 0.5 * sphericalbessely(2, k0_r)), ...
                                        struct('Index', 1, 'SiteOperator', a', 'SiteLabel', 1), ...
                                        struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
                                        correlations_Left, 2, 'unitary'});
%     int_type{end+1} = containers.Map(keySet, {Gamma_0 * (sphericalbesselj(0, k0_r) - 0.5 * sphericalbesselj(2, k0_r)), ...
%                                         struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
%                                         struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
%                                         correlations_Right, 1, 'dissipative_ARB'});
% 	int_type{end+1} = containers.Map(keySet, {Gamma_0 * (sphericalbesselj(0, k0_r) - 0.5 * sphericalbesselj(2, k0_r)), ...
%                                         struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
%                                         struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
%                                         correlations_Right, 1, 'dissipative_BRA'});
%     int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (sphericalbessely(0, k0_r) - 0.5 * sphericalbessely(2, k0_r)), ...
%                                         struct('Index', 1, 'SiteOperator', a, 'SiteLabel', clustersize), ...
%                                         struct('Index', 1, 'SiteOperator', a', 'SiteLabel', 1), ...
%                                         correlations_Right, 1, 'unitary'});
%     int_type{end+1} = containers.Map(keySet, {0.5 * Gamma_0 * (sphericalbessely(0, k0_r) - 0.5 * sphericalbessely(2, k0_r)), ...
%                                         struct('Index', 1, 'SiteOperator', a', 'SiteLabel', clustersize), ...
%                                         struct('Index', 1, 'SiteOperator', a, 'SiteLabel', 1), ...
%                                         correlations_Right, 1, 'unitary'});
end
input.subinput{1}.interactions = CreateInteractions(input, int_type);

%% run simulation
result = TimeIter(input, rho);

% %% output results
% time_array = (1:Nt)*dt;
% for time = 1:Nt
%     sigmaxfinal1(:, time) = result.hist{1}{time}.sigmax;
% end
% figure(1)
% plot(time_array, sigmaxfinal1');

sigma(rcount, Deltacount) = imag(E * d_eg * result.hist{1}{end}.rho(2,1));

end

% hold all;
% Delta_array = (1:nDelta)*deltaDelta + Delta_min;
% plot(Delta_array, sigma(rcount, :))

end

%% plot the cross-section as a function of detuning and lattice spacing
figure(1)
Delta_array = (1:nDelta)*deltaDelta + Delta_min;
r_array = (1:nr)*deltar + r_min;
[X, Y] = meshgrid(Delta_array, r_array);
surf(X, Y, sigma);

%% fit the cross-section to a Lorentzian to get half-decay rates and lineshifts for each set of lattice spacing
for i = 1:nr
  outputdata{i} = lorentzfitMODIFIED(Delta_array, sigma(i,:), [], [], '3');
  decayrate(i) = 2.0*sqrt(outputdata{i}{2}(3));
  lineshift(i) = outputdata{i}{2}(2);
end
% plot results
figure(2)
plot(r_array, decayrate, r_array, lineshift)
legend('half-decay rate', 'line shift')
xlabel('a / \lambda')
ylabel('Units of \Gamma_0')

toc