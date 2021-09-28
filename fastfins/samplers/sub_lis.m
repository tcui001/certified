function out = sub_lis(redu_def, model, prior, obs)
%SUB_LIS
%
% Builds likelihood informed subspace using subspace MCMC sampler
% DILI kernels are used
%
% Tiangang Cui, 27/Jan/2014
%
% The computational parameters y = RP\z, where inv(C) = RP'*RP,
% all the operation are defined on the computational parameters
%
% Note that when feed in the initial conditions, we need to transfrom
% the parameters onto y
%
% LIS:          the likelihood induced subspace
% CS:           the complement of the likelihood induced subspace
%
% Reference:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For tuning parameters of mcmc_def, see dili_mcmc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tuning parameters for compute the LIS, in the pass structure redu_def
%
% param_redu:
%       Used supplied reduced parameter space, if not provided, or set to
%       [], the algorithm will build one from scratch
%
% hess_type:
%       'PM' or 'GN', 'PM' for true Hessian with positive and negative
%       eigenvalues, 'GN' for Gauss Newton Hessian, default is 'GN'.
%
% gsvd_Nmax:
%       Max samples for building the LIS, default is 200
%
% gsvd_Nmin:
%       Min samples for building the LIS, default is 100
%
% gsvd_Ninter:
%       Number of steps for choosing the next sample to build LIS, default
%       is 200
%
% eigen_Nmax:
%       Max number of eigenvalues for the local decomposition, default 500
%
% gsvd_conv_tol:
%       Tolerance for stoping the adaptation of LIS, compared to the
%       initial distance, default 1E-4
%
% gsvd_trunc_tol:
%       Truncation threshold for the LIS,
%       if gsvd_trunc_tol > 0, truncate the SVD by S < gsvd_trunc_tol,
%       if gsvd_trunc_tol < 0, truncate the SVD by the accumulative energy
%       cumsum(S)/sum(S) < abs(gsvd_trunc_tol)
%
% eigen_tol:
%       Truncation threshold of the local eigendecomposition of the Hessian,
%       default 0.01.
%
% debug_flag:
%       With this option on, will also return the sample trace during the
%       adaptation, and the final global SVD.
%
% Pass-in function required:
%
% generalized_eigen = @(y, tol, Nmax) generalized_eigen(mesh, obs, ...
%                     prior, param, FEM, y, Nmax);
%
% [V, d] = generalized_eigen(y, tol, Nmax)
%
% input:
%       y,    parameters equipped with identity prior covariance
%       tol,  truncation of the local eigendecomposion
%       Nmax, maximum number of eigenvalues and eigenvectors want to compute
%
% output:
%       V, the eigenvectors
%       d, the eigenvalues, sorted in descending order
%
% Let H denote the Hessian of the llkd, and C the prior covariance, then
% this function compute the eigendecomposition of C^{0.5} H C^{0.5}.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mcmc_def, redu_def]    = process_input(redu_def);

num_basis               = zeros(redu_def.gsvd_Nmax, 1);
param_dist              = zeros(redu_def.gsvd_Nmax, 1);
samples                 = zeros(redu_def.gsvd_Nmax, redu_def.np);
stop_iter               = redu_def.gsvd_Nmax;
new_sample              = redu_def.init;
new_sigma               = mcmc_def.sigma;

[gsvd,param_redu]       = redu_init(redu_def, new_sample); % at the MAP

for i = 1:redu_def.gsvd_Nmax
    
    sub_def             = mcmc_def;
    sub_def.cov_def     = 'Restart';
    sub_def.mean        = param_redu.P'*redu_def.init;
    sub_def.sample_num  = param_redu.DoF*2;
    sub_def.cov         = diag((param_redu.S.^2 + 1).^(-1));
    sub_def.init        = param_redu.P'*new_sample;    
    sub_def.sigma       = new_sigma;
    sub_def.nstep       = redu_def.gsvd_Ninter;
    
    % determine the gsvd_Ninter, for adaptive jumpsize
    %if i <= redu_def.dyna_N
    %    mcmc_def.nstep  = redu_def.gsvd_Ninter*(redu_def.dyna_N+1-i);
    %else
    %    mcmc_def.nstep  = redu_def.gsvd_Ninter;
    %end
    
    % redefine the problem
    reduced             = basis_LIS(prior, param_redu.P);
    sub_def.minus_log_post  = @(v) minus_log_post(model, obs, reduced, v); 
        
    sub_chain           = amcmc(sub_def); % run subchain
    
    param_redu_p        = param_redu;
    new_sample          = param_redu_p.P*sub_chain.v_samples(end,:)';
    new_sigma           = sub_chain.sigma(end);

    [gsvd, param_redu]  = redu_reinit(redu_def, gsvd, new_sample);
    
    % check for convergence of the global basis, if not update global basis
    d                   = dist_fm(param_redu_p, param_redu);
    param_dist(i)       = d;
    num_basis (i)       = size(param_redu.P,2);
    samples   (i, :)    = new_sample;
    
    ni = 10;
    if i < ni
        md  = sum(param_dist(1:i))/i;
    else
        md  = sum(param_dist(i-ni+1:i))/ni;
    end
    
    % recompute the tolerance
    if i == ni
        gsvd_conv_tol   = md * redu_def.gsvd_conv_tol;
        fprintf('Convergence tol: %E \n', gsvd_conv_tol);
    end
    
    fprintf('%5d%5d%5d\t%E\t%E \n', [gsvd.Nsample, length(gsvd.S), param_redu.DoF, d, md]);
    
    if md < redu_def.gsvd_conv_tol && gsvd.Nsample > redu_def.gsvd_Nmin
        stop_iter       = i;
        break;
    end
    
    if gsvd.Nsample > redu_def.gsvd_Nmax
        stop_iter       = i;
        break;
    end
    
    % initialized the chain for next iteration
    
end

out.param_redu  = param_redu;
out.param_dist  = param_dist;
out.num_basis   = num_basis;
out.stop_iter   = stop_iter;
out.samples     = samples;
out.sigma       = sub_chain.sigma(end);

if redu_def.debug_flag 
    out.gsvd    = gsvd;
end

disp('done adaptation')

out.sub_cov.mu  = zeros(param_redu.DoF,1);
out.sub_cov.C   = diag((param_redu.S.^2 + 1).^(-1));
out.sub_cov.num = stop_iter;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mcmc_def, redu_def] = process_input(redu_def)

redu_def                    = lis_input(redu_def);
%number of sample intervals for global svd
if ~isfield(redu_def, 'gsvd_Ninter')
    redu_def.gsvd_Ninter    = 200;
end
% initial burning steps
%if ~isfield(redu_def, 'dyna_N')
%    redu_def.dyna_N         = 5;
%end

mcmc_def.init               = redu_def.init;
mcmc_def.nstep              = redu_def.gsvd_Ninter;
mcmc_def.proposal           = 'MALA';
user.high_d                 = false;
user.rate                   = 0.58;
mcmc_def                    = mcmc_input(mcmc_def, user);

end
