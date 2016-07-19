%%% an m file to calculate the location of the zeros of the mean-field
%%% Louivillian for the interacting Bose Hubbard model.  This model
%%% exhibits bistability so there should be a degeneracy of the zero value
%%% eigenvalue here.  Where does this eigenvalue come from and how does it
%%% evolve?

clear all;
hold all;

tic
addpath('../dev')
addpath('../test')

% create a set of inputs
clustersize = 1;
onsitedim = 15;
J = 3.0;
Delta = 1.5;
U = 2.0;
Omega = 3.79;        % the transition away from bistability occurs at 3.79-3.8
gamma = 0.3;

% reset the system
varargin = {'clustersize', clustersize, ...
            'onsitedim', onsitedim, ...
            'J', J, ...
            'Delta', Delta, ...
            'U', U, ...
            'Omega', Omega, ...
            'gamma', gamma, ...
            'L', { @L0, @LMF } };
        
absTol = 1e-7;
iterTol = 1e-3;     % tolerance for iterations
input.exists = true;
input = SetupSystem(input, varargin);
input = SetupEnvironment(input, {});
input = SetupTimeIter(input, {'probelist', { @SaveDensity }});
rho = InitializeRho(input);
 rho{1}(1,1)=0;
 rho{1}(3,3)=1;
% rho(onsitedim,onsitedim)=1;

% do the iteration
input.Nt = round(20.0 / (input.gamma * input.dt));
result = TimeIter(input, rho);

% check that the values for the mean-field bosonic coherence s
% and occupations are consistent with the theoretical values
% for each site
a_loc = annihilation(onsitedim);
for i = 1:clustersize
    a_site = InsertOperator(a_loc, i, onsitedim, clustersize);

    % calculate the mean-field bosonic coherence
    a_mf = trace(a_site * result.rho{1});
    % and input it into the theoretical equation... it should
    % be zero
    c = -2.0 * (-1.0 * Delta + 0.5i * gamma) / U;
    % effective coupling to one site is multiplied by the
    % system coordination number
    Jeff = J * input.coordination;  
    % in the paper, the driving frequency term is not divided by 2
    F = 0.5 * Omega;
    amf_eq = @(a) ((F - Jeff * a) / (-1.0 * Delta + 0.5i * gamma) ...
                    * hypergeom([], [1+c, c'], 8.0 * abs((F - Jeff * a) / U)^2) ...
                    / hypergeom([], [c, c'], 8.0 * abs((F - Jeff * a) / U)^2)) - a;
    amf_eq(a_mf)

    % and the expected occupation (different from paper in the
    % hypergeometric function's third argument, but correct...)
    a_occ = abs(2.0 * (F - Jeff * a_mf) / U)^2 * ...
                gammacomplex(c) * gammacomplex(c') / (gammacomplex(c+1) * gammacomplex(c'+1)) * ...
                hypergeom([], [c+1, c'+1], 8.0 * abs((F - Jeff * a_mf) / U)^2) ...
                / hypergeom([], [c, c'], 8.0 * abs((F - Jeff * a_mf) / U)^2);
    if(abs(trace(a_site' * a_site * result.rho{1}) - a_occ) > iterTol || amf_eq(a_mf) > iterTol)
        'Error!'
    end
    
    % find the eigenvalues of the superoperator matrix
    eigvals = eig(CreateSuperoperatorMatrix(@L0, input, result) ...
                    + CreateSuperoperatorMatrix(@LMF, input, result));
    scatter(real(eigvals),imag(eigvals));
end

for time = 1:numel(result.hist{1})
    nfinal(:, time) = result.hist{1}{time}.n;
end
figure(2)
plot(nfinal');

toc