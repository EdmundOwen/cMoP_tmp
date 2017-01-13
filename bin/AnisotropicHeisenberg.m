%%% a script to simulate the anisotropic Heisenberg model in two dimensions

close all; clear all;

tic
addpath('../dev')

%% define input variables
clustersizex = 2;
clustersizey = 2;
clustersize = clustersizex * clustersizey;
onsitedim = 2;

Omega = 0.0;
gamma = 1.0;
J = [10.0, 0.9, 1.0] / 2.0;

%% define the annihilation, creation and Pauli matrices
a = annihilation(onsitedim);
sigmax = a + a'; 
sigmay = 1i * (a - a');
sigmaz = a * a' - a' * a;
sigma = {sigmax, sigmay, sigmaz};

%% setup the inputs
input.exists = true;
input = SetupSystem(input, {'clustersize', clustersize, 'onsitedim', onsitedim, ...
                            'gamma', gamma, 'L', { @L0, @LMF }});%, @LMF}});%, @LBTSS }});
input = SetupEnvironment(input, {});

%% create an initial density matrix and solution
rho{1} = zeros(input.M);
rho{1}(1,1) = 1.0;
rho{1} = InsertOperator((eye(2)+a')/sqrt(2), 1, onsitedim, clustersize) * rho{1} * InsertOperator((eye(2)+a)/sqrt(2), 1, onsitedim, clustersize);
rho{1} = InsertOperator((eye(2)+a')/sqrt(2), 4, onsitedim, clustersize) * rho{1} * InsertOperator((eye(2)+a)/sqrt(2), 4, onsitedim, clustersize);
solution = {};
solution.rho{1} = rho;

%% create a new Hamiltonian for the cluster
H0 = sparse(input.M, input.M);
% first, the onsite driving term
Hloc = 0.5 * Omega * sigma{1};
H0 = H0 + TensorAddOperator(Hloc, onsitedim, clustersize);
% and secondly, the intersite coupling through the anisotropic Heisenberg
% interaction 
Hintx = sparse(input.M, input.M);
Hinty = sparse(input.M, input.M);
% between neighbouring sites in the x direction
Hintxloc = 0.5 * (J(1) * kron(sigma{1}, sigma{1}) ...
                + J(2) * kron(sigma{2}, sigma{2}) ...
                + J(3) * kron(sigma{3}, sigma{3}));
Hintxrow = TensorAddOperator(Hintxloc, onsitedim, clustersizex-1);
for i = 1:clustersizey
    if clustersizex ~= 1
        Hintx = Hintx + kron(speye(onsitedim^(clustersizex*(i-1))), ...
                                kron(Hintxrow, speye(onsitedim^(clustersizex*(clustersizey-i)))));
    end
end
% and y
Hintyloc = 0.5 * (J(1) * kron(sigma{1}, kron(speye(onsitedim^(clustersizex-1)), sigma{1})) ...
                + J(2) * kron(sigma{2}, kron(speye(onsitedim^(clustersizex-1)), sigma{2})) ...
                + J(3) * kron(sigma{3}, kron(speye(onsitedim^(clustersizex-1)), sigma{3})));
Hintyrow = TensorAddOperator(Hintyloc, onsitedim, clustersizex);
for i = 1:clustersizey-1
    if clustersizey ~= 1
        Hinty = Hinty + kron(speye(onsitedim^(clustersizex*(i-1))), ...
                                kron(Hintyrow, speye(onsitedim^(clustersizex*(clustersizey-i-1)))));
    end
end
H0 = H0 + Hintx + Hinty;
% and insert the Hamiltonian into the input struct
input.subinput{1}.H0 = H0;

%% create a new set of anisotropic Heisenberg interactions
int_list = {};
% the correlations between the interactions
correlationx = kron(ones(clustersizex), eye(2));
correlationy = kron(ones(clustersizey), eye(2));
correlation = [correlationx, zeros(2*clustersizex, 2*clustersizey);...
               zeros(2*clustersizey, 2*clustersizex), correlationy];
correlation = kron(ones(3), correlation);
           
% Pauli x, y, then z operators
count = 1;
for  i = 1:3
    % for interactions in the x direction
    for j = 1:clustersizey
        % to the left
        int_list{count} = {0.5 * J(i), (j-1)*clustersizex + 1, sigma{i}, ...
                            j*clustersizex, sigma{i}, correlation(count, :)};
        count = count + 1;
        % to the right
        int_list{count} = {0.5 * J(i), j*clustersizex, sigma{i}, ...
                            (j-1)*clustersizex + 1, sigma{i}, correlation(count, :)};
        count = count + 1;
    end
    % and for interactions in the y direction
    for j = 1:clustersizex
        % up
        int_list{count} = {0.5 * J(i), j, sigma{i}, ...
                            (clustersize-clustersizex) + j, sigma{i}, correlation(count, :)};
        count = count + 1;
        % down
        int_list{count} = {0.5 * J(i), (clustersize-clustersizex) + j, sigma{i}, ...
                            j, sigma{i}, correlation(count, :)};
        count = count + 1;
    end
end
% create the interactions
input.interactions = CreateInteractions(input, int_list);

% couple the partitions (redundant here)
for j = 1:numel(input.interactions)
	input.interactions{j}.B.Index = 1;
end
input.subinput{1}.interactions = input.interactions;

%% time iterate the original density matrix
time_varargin = {'Nt', 100, 'dt', 0.01, 'probelist', { @SaveDensity, @SaveCorrelation, @SaveMagnetisation }};
input = SetupTimeIter(input, time_varargin);
results = TimeIter(input, rho);

%% analyse the data ...
for time = 1:numel(results.hist{1})
    sigmaxfinal(:, time) = results.hist{1}{time}.sigmax;
end
figure(1)
plot(sigmaxfinal');
ylim([-1.0 1.0]);

for time = 1:numel(results.hist{1})
    sigmayfinal(:, time) = results.hist{1}{time}.sigmay;
end
figure(2)
plot(sigmayfinal');
ylim([-1.0 1.0]);

for time = 1:numel(results.hist{1})
    sigmazfinal(:, time) = results.hist{1}{time}.sigmaz;
end
nsteady = results.hist{1}{time}.n;
figure(3)
plot(sigmazfinal');
ylim([-1.0 1.0]);

toc