close all; clear all;

addpath('../dev')

tic
% for m = 21:30
%     for n = 1:10
gamma = 0.0;
Omega = 0.0;
setup_varargin = {'clustersize', 7, 'onsitedim', 2, 'L', { @L0 }, 'gamma', gamma, 'Omega', Omega};
timeiter_varargin = { 'method', 'heun', 'Nt', 1000, 'dt', 0.01, 'probelist', { @SaveDensity, @SaveCorrelation } };

input.exists = true;
input = SetupSystem(input, setup_varargin);
input = SetupTimeIter(input, timeiter_varargin);

onsitedim = input.onsitedim;
clustersize = input.clustersize;

a = annihilation(onsitedim);
%% edit this
wavefn = sparse(onsitedim^clustersize, 1);
wavefn(1) = 1.0;

a_wann1 = (kron(a, speye(onsitedim^(clustersize - 1))) ...
            - sqrt(2) * kron(kron(speye(onsitedim), a), speye(onsitedim^(clustersize - 2))) ...
            + kron(kron(speye(onsitedim^2), a), speye(onsitedim^(clustersize - 3)))) / 2;
        
a_wann2 = (kron(kron(speye(onsitedim^(2)), a), speye(onsitedim^4)) ...
            - sqrt(2) * kron(kron(speye(onsitedim^(3)), a), speye(onsitedim^3)) ...
            + kron(kron(speye(onsitedim^(4)), a), speye(onsitedim^2))) / 2;
        
a_wann3 = (kron(kron(speye(onsitedim^(4)), a),speye(onsitedim^2) )...
            - sqrt(2) * kron(kron(speye(onsitedim^(5)), a), speye(onsitedim)) ...
            + kron(kron(speye(onsitedim^(6)), a), speye(onsitedim^0))) / 2;
        
% a_wann2 = (kron(speye(onsitedim^(clustersize - 1)), a) ...
%             - sqrt(2) * kron(kron(speye(onsitedim^(clustersize-2)), a), speye(onsitedim)) ...
%             + kron(kron(speye(onsitedim^(clustersize-3)), a), speye(onsitedim^2))) / 2;

wavefn = (0.8*a_wann1' + 0.2 * a_wann2') * wavefn;
%wavefn = kron(a', speye(onsitedim^(clustersize-1))) * wavefn;
rho{1} = kron(wavefn', wavefn);

t = 1.0;
tprime = sqrt(2) * t;
delta = 2.0 * t;
U = 5.0 * t;


H = tprime * kron(kron(a, a'), speye(onsitedim^(clustersize - 2)));
for i = 2:clustersize-1
    H = H + tprime * kron(kron(kron(speye(onsitedim^(i - 1)), a), a'), speye(onsitedim^(clustersize - i - 1)));
end

for i = 1:(clustersize-1)/2-1
    H = H + t * kron(kron(kron(kron(speye(onsitedim^(i * 2 - 1)), a), ...
                speye(onsitedim)), a'), speye(onsitedim^(clustersize - 2 * (i + 1))));
end

% for i = 1:(clustersize-1)/2
%     H = H + t * kron(kron(kron(kron(speye(onsitedim^(i * 2 - 2)), a), ...
%                 speye(onsitedim)), a'), speye(onsitedim^(clustersize - 2 * (i + 1) + 1)));
% end

pattern = [0.5, -1.0/sqrt(2), 0.5, 0, 0.5, -1.0/sqrt(2), 0.5];

% factor of 0.5 for delta is due to H + H'
for i = 1:clustersize
    H = H + kron(kron(speye(onsitedim^(i-1)), 0.5 * delta * (a') * a ...
                                            + 0.5 * U * (a') * (a') * a * a ...
                                            + 0.5 * input.Omega * a), speye(onsitedim^(clustersize - i)));
end

input.H0 = H + H';

results = TimeIter(input, rho);

Nt = input.Nt;
ntmp = results.hist{1}{Nt}.n;

% dens(n,m) = sum(ntmp);
% dens12(n,m) = ntmp(1)-ntmp(2);
% dens23(n,m) = ntmp(2)-ntmp(3);
% dens24(n,m) = ntmp(2)-ntmp(4);
% dens34(n,m) = ntmp(3)-ntmp(4);
% 
% %     end
% % end
% 
% figure(1)
% pcolor(dens)
% figure(2)
% pcolor(dens12)
% figure(3)
% pcolor(dens23)
% figure(4)
% pcolor(dens24)
% figure(5)
% pcolor(dens34)

for time = 1:numel(results.hist{1})
	nfinal(:, time) = results.hist{1}{time}.n;
end
figure(1)
plot(nfinal')

for time = 1:numel(results.hist{1})
    g2final(:, time) = results.hist{1}{time}.g2;
end
figure(2)
plot(real(g2final'))
figure(3)
plot(imag(g2final'))

[vecs vals] = eig(full(results.rho{1}));
for index = 1:onsitedim^clustersize
    tmp.rho{1} = kron(vecs(:,index)', vecs(:,index));
    tmp = SaveDensity(tmp, input.subinput{1}, index);
end
for i = 1:onsitedim^clustersize
    n(:,i)=tmp.hist{1}{i}.n * vals(i,i);
end
figure(4)
plot(n)
hold
plot(nfinal(:,end), 'linewidth', 2);
hold off

toc