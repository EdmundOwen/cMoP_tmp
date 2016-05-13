%%% a script to calculate the steady state of a driven dissipative sawtooth
%%% lattice using a Wannier state representation of the flat band with the
%%% higher bands traced out

tic
addpath('../dev')

%% define input variables
clustersize = 2;
onsitedim = 2;
t = 1.0;
U = 2.0;
OmegaA = 0.01;
OmegaB = 0.01;
gamma = 0.001;

%% setup the inputs
input.exists = true;
input = SetupSystem(input, {'clustersize', clustersize, 'onsitedim', onsitedim, 'gamma', gamma, 'operators', { @L0, @LMF, @LBTSS }});
input = SetupEnvironment(input, {});
input = SetupSteadyState(input, {});

%% create an initial density matrix and solution
rho = zeros(input.M);
rho(1,1) = 1.0;
solution = {};
solution.rho = rho;

%% define the Wannier state coefficients
wA0 = 2.026 / (2.0 * pi)^1.5;
wA1 = -0.253 / (2.0 * pi)^1.5;
wB0 = -4.686 / (2.0 * pi)^1.5;
wB1 = 0.634 / (2.0 * pi)^1.5;

%% create a new Hamiltonian for the Wannier states
H0 = sparse(input.M, input.M);
a = annihilation(input.onsitedim);
% first, the onsite energy, onsite interaction and driving term
Hloc = 0.5 * U * (2.0 * wA0^4 + wB0^4) * (a')*(a')*a*a ...
        + 0.5 * (2.0 * OmegaA * wA0 + OmegaB * wB0) * (a') ...
        + 0.5 * (2.0 * OmegaA' * wA0 + OmegaB' * wB0) * a;
H0 = H0 + TensorAddOperator(Hloc, onsitedim, clustersize);
% and secondly, the intersite coupling through the onsite A interaction
Hint = 0.5 * U * wA0^4 * (kron((a * a)', a * a) + kron(a * a, (a * a)') ...
                         + 4.0 * kron((a') * a, (a') * a) ...
                         + 2.0 * (kron(a', (a') * a * a) + kron((a') * a * a, a')));
H0 = H0 + TensorAddOperator(Hint, onsitedim, clustersize-1);
% and insert the Hamiltonian into the input struct
input.H0 = H0;

%% create a new set of interactions for the Wannier states
int_list = {};
fac = 0.5 * U * wA0^4;
% the correlations between the interactions
correlation = kron(eye(2), ones(7));
int_list{1} = {fac, 1, a*a, clustersize, (a*a)', correlation(1, :)};
int_list{2} = {fac, 1, (a*a)', clustersize, a*a, correlation(2, :)};
int_list{3} = {4.0 * fac, 1, (a')*a, clustersize, (a')*a, correlation(3, :)};
int_list{4} = {2.0 * fac, 1, (a')*a*a, clustersize, a', correlation(4, :)};
int_list{5} = {2.0 * fac, 1, a', clustersize, (a')*a*a, correlation(5, :)};
int_list{6} = {2.0 * fac, 1, (a')*a*(a'), clustersize, a, correlation(6, :)};
int_list{7} = {2.0 * fac, 1, a, clustersize, (a')*a*(a'), correlation(7, :)};
int_list{8} = {fac, clustersize, a*a, 1, (a*a)', correlation(8, :)};
int_list{9} = {fac, clustersize, (a*a)', 1, a*a, correlation(9, :)};
int_list{10} = {4.0 * fac, clustersize, (a')*a, 1, (a')*a, correlation(10, :)};
int_list{11} = {2.0 * fac, clustersize, (a')*a*a, 1, a', correlation(11, :)};
int_list{12} = {2.0 * fac, clustersize, a', 1, (a')*a*a, correlation(12, :)};
int_list{13} = {2.0 * fac, clustersize, (a')*a*(a'), 1, a, correlation(13, :)};
int_list{14} = {2.0 * fac, clustersize, a, 1, (a')*a*(a'), correlation(14, :)};
% create the interactions
input = CreateInteractions(input, int_list);

%% find the steady state solution
input.verbose = true;
result = CalculateSteadyState(input, rho, solution);

%% analyse the data...
g2 = kron((a') * a, (a') * a);
n2 = kron(a' * a, eye(onsitedim));
nval = trace(n2 * result)
g2val = trace(g2 * result) / nval^2

toc