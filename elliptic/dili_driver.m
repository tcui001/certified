% This is a setup file for DILI options
% Tiangang Cui, 01/Oct/2015

mcmc_def.init        = vmap;    % MAP estimate as the initial guess
mcmc_def.nstep       = 1E5;     % total number of steps
mcmc_def.save_batch  = 1;       % save the results for every nbatch number of iterations
mcmc_def.sigma       = -1;      % subspace jump size
mcmc_def.using_gibbs = true;    % Gibbs update
mcmc_def.proposal    = 'MALA';  % Or MALA for MGLI-Prior proposal
mcmc_def.dt          = 2;       % complement space jump size
mcmc_def.projection  = [];      % a projection matrix for saving the subspace MCMC samples 

% posterior density function
mcmc_def.minus_log_post     = @(v) minus_log_post(model, obs, prior, v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define adaptation schemes
redu_def.init           = vmap; 
redu_def.eigen_Nmax     = obs.Ndata;    % maximumn number of eigenvector for each local eigendecomposition
redu_def.eigen_tol      = 1E-2;         % local trucation threshold
redu_def.gsvd_trunc_tol = 1E-1;          % global truncation threshold

redu_def.gsvd_Nmax      = 500;          % max number of samples used for computing the LIS
redu_def.gsvd_Nmin      = 100;          % min number of samples for computing the LIS
redu_def.gsvd_Ninter    = 100;          % lag for selecting samples
redu_def.gsvd_conv_tol  = 1E-5;         % relative convergence tol for the LIS
redu_def.method         = 'Eig';        % Using eigendecompostion rather than SVD for building the local LIS

% functions
redu_def.PPGNH          = @(v, tol, Nmax) eigen_PPGNH(model, obs, prior, v, tol, Nmax);

tic;
dimredu             = dili_lis(mcmc_def, redu_def);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


param_redu          = dimredu.param_redu;

mcmc_def.save_all_flag      = true;    % save all the parameters
mcmc_def.projection         = param_redu.P;      % projection matrix for the marginal MCMC history
mcmc_def.nstep              = 1E6;

tic
out_dili                    = dili_mcmc(mcmc_def, param_redu);
toc