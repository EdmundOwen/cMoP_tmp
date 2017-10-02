%%% a script to simulate the anisotropic Heisenberg model in two dimensions

close all; clear all;

tic
addpath('../dev')

%% define input variables
clustersizex = 1;
clustersizey = 1;
clustersize = clustersizex * clustersizey;
onsitedim = 2;
latticedim = 100;

gamma = 1.0;
Delta = 0.0;
J = [-7.0, 6.0, 2.0]; Omega = 1.0;

%% define the annihilation, creation and Pauli matrices
a = annihilation(onsitedim);
sigmax = a + a'; 
sigmay = -1i * (a - a');
sigmaz = a * a' - a' * a;
sigma = {sigmax, sigmay, sigmaz};

%% setup the inputs
input.exists = true;
input = SetupSystem(input, {'clustersize', clustersize, 'onsitedim', onsitedim, ...
                            'gamma', gamma, 'Omega', Omega, 'Delta', Delta, 'L', { @L0, @LMF }});%, @LMF}});%, @LBTSS }});
                        
%% HACK!!!  The raising spin operator is the transpose of the QHO creation operator (and similarly for the lowering operator)
input.A_Lindblad{1} = [0, 0; 1, 0];
                        
input = SetupEnvironment(input, {});
% set interactions to null as this is a partitioned system and this
% shouldn't be used
input.interactions = [];

%% setup partitions
input.noPartitions = 2;
input.subinput{1} = SetupPartition(input, 1);
input.subinput{2} = SetupPartition(input, 2);

%% create an initial density matrix and solution
rhovec1 = [0.3575, 0.5155, 0.4863];
rhovec2 = [-0.2155, 0.3110, -0.6576];

rho{1} = 0.5*(eye(2)+rhovec1(1)*sigma{1}+rhovec1(2)*sigma{2}+rhovec1(3)*sigma{3});
rho{2} = 0.5*(eye(2)+rhovec2(1)*sigma{1}+rhovec2(2)*sigma{2}+rhovec2(3)*sigma{3});

solution = {};
solution.rho = rho;

%% create a new set of anisotropic Heisenberg interactions
int_list_1 = cell(1,2);
int_list_2 = cell(1,2);

keySet = {'interactionStrength', 'A', 'B', 'correlations', 'coordination'};
% Pauli x, y, then z operators
for  i = 1:3
    % for the first partition
    int_list_1{i} = containers.Map(keySet, ...
                    {0.5 * J(i) / latticedim, struct('Index', 1, 'SiteOperator', sigma{i}, 'SiteLabel', 1), ...
                         struct('Index', 2, 'SiteOperator', sigma{i}, 'SiteLabel', 1), [1 1 1], 2 * latticedim});
    % and for the second partition
    int_list_2{i} = containers.Map(keySet, ...
                    {0.5 * J(i) / latticedim, struct('Index', 2, 'SiteOperator', sigma{i}, 'SiteLabel', 1), ...
                         struct('Index', 1, 'SiteOperator', sigma{i}, 'SiteLabel', 1), [1 1 1], 2 * latticedim});
end
input.subinput{1}.interactions = CreateInteractions(input, int_list_1);
input.subinput{2}.interactions = CreateInteractions(input, int_list_2);

% tic
% input = SetupSteadyState(input, {'SSError', 1e-5, 'Niter', 5000});
% input.verbose = false;
% result = CalculateSteadyState(input, rho, solution);
% toc 
% return

%% time iterate the original density matrix
Nt = 5001; dt = 0.01; time_array = (1:Nt) * dt;
time_varargin = {'Nt', Nt, 'dt', dt, 'probelist', { @SaveRho, @SaveDensity, @SaveMagnetisation }};
input = SetupTimeIter(input, time_varargin);
results = TimeIter(input, rho);

%% analyse the data ...
for time = 1:Nt
    sigmaxfinal1(:, time) = results.hist{1}{time}.sigmax;
    sigmaxfinal2(:, time) = results.hist{2}{time}.sigmax;
end
figure(1)
plot(time_array, sigmaxfinal1', time_array, sigmaxfinal2');
ylim([-1.0 1.0]);

for time = 1:Nt
    sigmayfinal1(:, time) = results.hist{1}{time}.sigmay;
    sigmayfinal2(:, time) = results.hist{2}{time}.sigmay;
end
figure(2)
plot(time_array, sigmayfinal1', time_array, sigmayfinal2');
ylim([-1.0 1.0]);

for time = 1:Nt
    sigmazfinal1(:, time) = results.hist{1}{time}.sigmaz;
    sigmazfinal2(:, time) = results.hist{2}{time}.sigmaz;
end
figure(3)
plot(time_array, sigmazfinal1', time_array, sigmazfinal2');
ylim([-1.0 1.0]);

toc