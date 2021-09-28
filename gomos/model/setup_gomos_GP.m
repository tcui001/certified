function [model, obs, prior] = setup_gomos_GP()

% this script setup the GOMOS example with Gaussian process prior
%
% Tiangang Cui, 18/May/2014

%**************************************************************************
% load data sets
%**************************************************************************

load('transmission50.mat');
% Loads: alts nalts wavel nwavel ngas C A X T weights
%
% alts - altitudes [km]
% wavel - wavelengts [nm]
% C - cross sections, nwavel*ngas
% A - geometry matrix, nalts*nalts
% X - gas profiles, nalts*ngas (ngas*nalts in inversion)
% T - true transmissions, nwavel*nalts
% weights - knonw std error in transmissions, nwavel*ngas
%
% Direct model: T = exp(-C*X'*A')
% Observations Tobs = T + randn(size(T))*weights

% generate observations
% Tobs = T + randn(size(T)).*weights;

load('Tobs.mat')
% Scale C for numerical stability
scale           = mean(C); 
C               = bsxfun(@rdivide,C,scale);

%  scaled true profiles
Xs              = bsxfun(@times,X,scale);

% model function
% modelfun      = @(x) exp(-C*reshape(x,ngas,nalts)*A');
AC              = kron(A,C); % needed for Jacobian calculations, (nalts*nwavel)*(nalts*ngas)

% prior
sig             = (1*std(log(Xs))).^2;               % variances
L               = (10*ones(1,4)).^2;                 % correlation lengths ^2
tau             = 1e-9*ones(1,4);                    % nugget
%
mu              = mean(log(Xs));                       % means
mu              = (repmat(mu,nalts,1))';
mu              = mu(:);

%**************************************************************************
% setup prior, model and observation data structure
%**************************************************************************

Gamma           = gomos_prior_GP(alts,sig,L,tau);    % prior cov. matrix
prior           = set_prior_GP(Gamma, mu);           % wrap around the prior
prior.true_x    = Xs;

model.A         = A;
model.C         = C;
model.AC        = AC;
model.nalts     = nalts;
model.ngas      = ngas;
model.nwave     = length(wavel);
model.wavel     = wavel;
model.T         = T;
model.alts      = alts;
model.scale     = scale;

obs.std         = weights(:);
obs.data        = Tobs(:);
obs.Ndata       = size(obs.data, 1);
obs.Nsensors    = obs.Ndata;
obs.Ndatasets   = 1;

end