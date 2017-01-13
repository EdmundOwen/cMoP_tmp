%%% a script to simulate the anisotropic Heisenberg model in two dimensions

close all; clear all;

tic
addpath('../dev')

%% define input variables
clustersizex = 1;
clustersizey = 1;
clustersize = clustersizex * clustersizey;
onsitedim = 2;

Omega = 1.0;
gamma = 1.0;
%J = [-5.5, 6.0, 2.0] * 2.0;  % multiply by 2 to get the coordination number
J = [-7.0, 6.0, 2.0];  % multiply by 2 to get the coordination number

%% define the annihilation, creation and Pauli matrices
a = annihilation(onsitedim);
sigmax = a + a'; 
sigmay = 1i * (a - a');
sigmaz = a * a' - a' * a;
sigma = {sigmax, sigmay, sigmaz};

%% setup the inputs
input.exists = true;
input = SetupSystem(input, {'clustersize', clustersize, 'onsitedim', onsitedim, ...
                            'gamma', gamma, 'Omega', Omega, 'Delta', 0.0, 'L', { @L0, @LMF }});%, @LMF}});%, @LBTSS }});
input = SetupEnvironment(input, {});

%% setup partitions
input.noPartitions = 2;
input.subinput{1} = SetupPartition(input, 1);
input.subinput{2} = SetupPartition(input, 2);

%% create an initial density matrix and solution
rho{1} = zeros(input.M);
rho{1}(1,1) = 1.0;
rho{2} = InsertOperator((eye(2)+a')/sqrt(2), 1, onsitedim, clustersize) * rho{1} * InsertOperator((eye(2)+a)/sqrt(2), 1, onsitedim, clustersize);
solution = {};
solution.rho = rho;

%% create a new set of anisotropic Heisenberg interactions
int_list = {};
% set up the couplings between the partitions
input.subinput{1}.interactions = {};
input.subinput{2}.interactions = {};

% Pauli x, y, then z operators
count = 1;
for  i = 1:3
    % for the first partition
    input.subinput{1}.interactions{end+1} ...
         = struct('interactionStrength', 0.5 * J(i), 'A', struct('Operator', sigma{i}), ...
                    'B', struct('Index', 2, 'Operator', sigma{i}), 'Correlations', [1 1 1]);
    % and for the second partition
    input.subinput{2}.interactions{end+1} ...
         = struct('interactionStrength', 0.5 * J(i), 'A', struct('Operator', sigma{i}), ...
                    'B', struct('Index', 1, 'Operator', sigma{i}), 'Correlations', [1 1 1]);
end

%% time iterate the original density matrix
time_varargin = {'Nt', 10000, 'dt', 0.005, 'probelist', { @SaveDensity, @SaveCorrelation, @SaveMagnetisation }};
input = SetupTimeIter(input, time_varargin);
results = TimeIter(input, rho);

%% analyse the data ...
for time = 1:numel(results.hist{1})
    sigmaxfinal1(:, time) = results.hist{1}{time}.sigmax;
    sigmaxfinal2(:, time) = results.hist{2}{time}.sigmax;
end
figure(1)
plot(1:time, sigmaxfinal1', 1:time, sigmaxfinal2');
ylim([-1.0 1.0]);

for time = 1:numel(results.hist{1})
    sigmayfinal1(:, time) = results.hist{1}{time}.sigmay;
    sigmayfinal2(:, time) = results.hist{2}{time}.sigmay;
end
figure(2)
plot(1:time, sigmayfinal1', 1:time, sigmayfinal2');
ylim([-1.0 1.0]);

for time = 1:numel(results.hist{1})
    sigmazfinal1(:, time) = results.hist{1}{time}.sigmaz;
    sigmazfinal2(:, time) = results.hist{2}{time}.sigmaz;
end
figure(3)
plot(1:time, sigmazfinal1', 1:time, sigmazfinal2');
ylim([-1.0 1.0]);

toc