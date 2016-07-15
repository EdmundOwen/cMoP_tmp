%%% a script to calculate the steady state of a driven dissipative sawtooth
%%% lattice using a Wannier state representation of the flat band with the
%%% higher bands traced out

close all; clear all;

tic
addpath('../dev')
% 
% nOmega=5;
% dOmega=1.0;
% 
% for omegastep = 0:nOmega
%% define input variables
clustersize = 4;
onsitedim = 2;
t = 1.0;
U = 0.0;
OmegaA = 0.1;
OmegaB = 1.0;

% Omega = omegastep*dOmega;
Omega = OmegaA;
delta = 0.1;
gamma = 1.0;

gammaA = gamma;
gammaB = 0.0;
OmegaEnd = delta;

%% setup the inputs
input.exists = true;
input = SetupSystem(input, {'clustersize', clustersize, 'onsitedim', onsitedim, ...
                            'gamma', gamma, 'operators', { @L0 }});%, @LMF}});%, @LBTSS }});
input = SetupEnvironment(input, {});
input = SetupSteadyState(input, {});

%% create an initial density matrix and solution
rho = num2cell(zeros(input.M));
rho{1}(1,1) = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assym = 0.3;
a = annihilation(input.onsitedim);
wavefn = zeros(input.M, 1);
wavefn(1) = 1.0;
% wavefn = InsertOperator(a', 1, onsitedim, clustersize) * InsertOperator(a', 3, onsitedim, clustersize) * wavefn;
rho{1} = kron(wavefn', wavefn);
% rho = (0.5 + assym) * rho + (0.5 - assym) * InsertOperator(a', 1, onsitedim, clustersize) * rho * InsertOperator(a, 1, onsitedim, clustersize);
% rho = (0.5 - assym) * rho + (0.5 + assym) * InsertOperator(a', 2, onsitedim, clustersize) * rho * InsertOperator(a, 2, onsitedim, clustersize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solution = {};
solution.rho{1} = rho;

%% define the Wannier state coefficients
wA0 = 0.456018;%2.865 / (2.0 * pi)^1.;
wA1 = -0.0569315;%-0.253 / (2.0 * pi)^1.;
wB0 = -0.745749;%-4.686 / (2.0 * pi)^1.;
wB1 = 0.100842;%0.634 / (2.0 * pi)^1.;

%% renormalise the input parameters to remove the wA0 and wB0 dependence
gammaA = gammaA / wA0^2;
gammaB = gammaB / wB0^2;
OmegaA = 0.0;
OmegaB = Omega / wB0;

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
% add an additional driving term on the first site
H0 = H0 + InsertOperator(0.5 * (OmegaEnd * (a') + OmegaEnd' * a), 1, onsitedim, clustersize);
% and insert the Hamiltonian into the input struct
input.subinput{1}.H0 = H0;

%% alter the Lindblad weights to make them applicable for the system Wannier states
% diagonal contribution
input.subinput{1}.Lindblad_weights = (2.0 * gammaA * wA0^2 + gammaB * wB0^2) * eye(clustersize);
% off-diagonal contributions
for i = 1:clustersize-1
    input.subinput{1}.Lindblad_weights(i, i+1) = gammaA * wA0^2;
    input.subinput{1}.Lindblad_weights(i+1, i) = input.Lindblad_weights(i, i+1);
end
input.subinput{1}.Lindblad_weights(1,end)=input.Lindblad_weights(1,2);
input.subinput{1}.Lindblad_weights(end,1)=input.Lindblad_weights(end,end-1);

%% create a new set of interactions for the Wannier states
int_list = {};
fac = 0.5 * U * wA0^4;
% the correlations between the interactions
correlation = kron(eye(2), ones(8));
int_list{1} = {fac, 1, a*a, clustersize, (a*a)', correlation(1, :)};
int_list{2} = {fac, 1, (a*a)', clustersize, a*a, correlation(2, :)};
int_list{3} = {4.0 * fac, 1, (a')*a, clustersize, (a')*a, correlation(3, :)};
int_list{4} = {2.0 * fac, 1, (a')*a*a, clustersize, a', correlation(4, :)};
int_list{5} = {2.0 * fac, 1, a', clustersize, (a')*a*a, correlation(5, :)};
int_list{6} = {2.0 * fac, 1, (a')*a*(a'), clustersize, a, correlation(6, :)};
int_list{7} = {2.0 * fac, 1, a, clustersize, (a')*a*(a'), correlation(7, :)};
int_list{8} = {gamma * wA0^2, 1, a, clustersize, a, correlation(8,:), 'dissipative'};
int_list{9} = {fac, clustersize, a*a, 1, (a*a)', correlation(9, :)};
int_list{10} = {fac, clustersize, (a*a)', 1, a*a, correlation(10, :)};
int_list{11} = {4.0 * fac, clustersize, (a')*a, 1, (a')*a, correlation(11, :)};
int_list{12} = {2.0 * fac, clustersize, (a')*a*a, 1, a', correlation(12, :)};
int_list{13} = {2.0 * fac, clustersize, a', 1, (a')*a*a, correlation(13, :)};
int_list{14} = {2.0 * fac, clustersize, (a')*a*(a'), 1, a, correlation(14, :)};
int_list{15} = {2.0 * fac, clustersize, a, 1, (a')*a*(a'), correlation(15, :)};
int_list{16} = {gamma * wA0^2, clustersize, a, 1, a, correlation(16,:), 'dissipative'};
% create the interactions
interactions = CreateInteractions(input, int_list);

%% Add the set of dissipator interactions.  These are due to the non-local
% nature of the dissipation in the Wannier basis
% interactions{8} = {0.5i * gamma * wA0^2, input.A_Lindblad{2}, input.A_Lindblad{1}', correlation(8,:)};
% interactions{9} = {-0.5i * gamma * wA0^2, input.A_Lindblad{2}', input.A_Lindblad{1}, correlation(9,:)};
% interactions{17} = {0.5i * gamma * wA0^2, input.A_Lindblad{1}, input.A_Lindblad{2}', correlation(17,:)};
% interactions{18} = {-0.5i * gamma * wA0^2, input.A_Lindblad{1}', input.A_Lindblad{2}, correlation(18,:)};
% set the interactions in the input struct
input.subinput{1}.interactions = interactions;

% %% find the steady state solution
% input.verbose = true;
% result = CalculateSteadyState(input, rho, solution);
% 
% %% analyse the data...
% g2 = kron(kron((a') * a, (a') * a), eye(onsitedim^(clustersize-2)));
% n2 = kron(a' * a, eye(onsitedim^(clustersize-1)));
% nval = trace(n2 * result)
% g2val = trace(g2 * result) / nval^2

%% time iterate the original density matrix
time_varargin = {'Nt', 1000, 'dt', 0.01, 'probes', { @SaveDensity, @SaveCorrelation }};
input = SetupTimeIter(input, time_varargin);
results = TimeIter(input, rho);

%% analyse the data ...
for time = 1:numel(results.hist{1})
    nfinal(:, time) = results.hist{1}{time}.n;
end
figure(1)
plot(nfinal');
% ylim([0 0.5]);

for time = 1:numel(results.hist{1})
    g2final(:, time) = results.hist{1}{time}.g2;
end

figure(2)
plot(real(g2final'))
figure(3)
plot(imag(g2final'))
% 
% fprintf('i=%d,  n=%d\n', omegastep, sum(nfinal(:,end))/clustersize)
% 
% dens{omegastep+1}=[omegastep*dOmega, sum(nfinal(:,end))/clustersize];
% end
% 
% for i = 0:nOmega
%     omega(i+1)=dens{i+1}(1);
%     densval(i+1)=dens{i+1}(2);
% end
% plot(omega,densval)
% save 'matlabdens_d2N6' dens

toc