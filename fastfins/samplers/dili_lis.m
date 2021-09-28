function out = dili_lis(mcmc_def, redu_def)
%DILI_LIS
%
% builds likelihood informed subspace using DILI MCMC sampler
%
% Tiangang Cui, 17/Jan/2014
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mcmc_def, redu_def, out] = process_input(mcmc_def, redu_def);

sigma           = mcmc_def.sigma;
num_basis       = zeros(redu_def.gsvd_Nmax, 1);
param_dist      = zeros(redu_def.gsvd_Nmax,1);

[gsvd, param_redu]      = redu_init(redu_def, redu_def.init); % at the MAP
[curr, stat, kernel]    = dili_init(mcmc_def, param_redu, sigma); % initialization

%%%%%%%%%%%%%%%%%%
% MCMC
acc_sub         = 0;
acc_null        = 0;
batch           = 0;
stop_iter       = redu_def.gsvd_Nmax;
ii              = 0;

for i = 1:mcmc_def.nstep
    r               = randn(mcmc_def.np,1); % random number
    % project random number
    r_sub           = param_redu.P'*r;
    r_null          = r - param_redu.P*r_sub;
    % propose in the subspace, and evaluate the acceptance rate
    [alpha, next]   = propose_sub(mcmc_def, param_redu, kernel, curr, r_sub);
    if  log(rand)   < alpha
        acc_sub     = acc_sub+1;
        curr        = next; % update
    end
    
    if ~redu_def.sub_chain
        % propose in the null, and evaluate the acceptance rate
        [alpha, next]   = propose_null(mcmc_def, param_redu, kernel, curr, r_null);
        if  log(rand)   < alpha
            acc_null    = acc_null+1;
            curr        = next; % update
        end
    end
    
    % update the emperical mean and covariance in the LIS
    stat.cross  = stat.cross + curr.v_sub(:)*curr.v_sub(:)';
    stat.sum    = stat.sum   + curr.v_sub;
    stat.num    = stat.num   + 1;
    
    % save
    if  mcmc_def.save_batch  == 1
        out.j                       = out.j+1;
        out.v_samples(out.j,:)      = curr.v;
        out.llkd(out.j)             = curr.mllkd;
        out.lpt (out.j)             = curr.mlpt;
    elseif mod(i, mcmc_def.save_batch)   == 0
        out.j                       = out.j+1;
        out.v_samples(out.j,:)      = curr.v;
        out.llkd(out.j)             = curr.mllkd;
        out.lpt (out.j)             = curr.mlpt;
    end
        
    % tune the transition kernel
    batch       =  batch + 1;
    if  batch   == mcmc_def.nbatch
        delta   =  min(0.1,sqrt(mcmc_def.nbatch/i));
        if (acc_sub/mcmc_def.nbatch) < mcmc_def.rate
            sigma   = sigma - delta;
        else
            sigma   = sigma + delta;
        end
        
        % save
        out.k               = out.k+1;
        out.sigma(out.k)    = sigma;
        out.acc(out.k,:)    = [acc_sub, acc_null]/mcmc_def.nbatch;
        
        % reset coounter
        batch       = 0;
        acc_sub     = 0;
        acc_null    = 0;
        % build new kernel
        stat        = build_cov(mcmc_def, stat);
        kernel      = build_kernel(mcmc_def, param_redu, stat, sigma);
    end
    
    % determine the gsvd_Ninter, for adaptive jumpsize
    if i <= redu_def.dyna_N
        Ninter  = redu_def.gsvd_Ninter*(redu_def.dyna_N+1-i);
    else
        Ninter  = redu_def.gsvd_Ninter;
    end
    
    
    % rebuild the subspace
    if  mod(i,Ninter) == 0
        % update lis, re initilaize the MCMC
        param_redu_p            = param_redu;
        [gsvd, param_redu]      = redu_reinit(redu_def, gsvd, curr.v);
        if redu_def.sub_chain
            [curr,stat,kernel]  = dili_reinit(mcmc_def, curr, param_redu, sigma, [], 0);
        else
            [curr,stat,kernel]  = dili_reinit(mcmc_def, curr, param_redu, sigma, out.v_samples(1:out.j,:), mcmc_def.save_batch);
        end
                
        % check for convergence of the global basis, if not update global basis
        ii              = ii + 1;
        d               = dist_fm(param_redu_p, param_redu);
        param_dist(ii)  = d;
        num_basis(ii)   = size(param_redu.P,2);
        
        ni = 10;
        if ii < ni
            md = sum(param_dist(1:ii))/ii;
        else
            md = sum(param_dist(ii-ni+1:ii))/ni;
        end
        
        % recompute the tolerance
        if ii == ni
            gsvd_conv_tol   = md * redu_def.gsvd_conv_tol;
            fprintf('Convergence tol: %E\n\n', gsvd_conv_tol);
        end
        
        fprintf('%5d%5d%5d\t%E\t%E\n', [gsvd.Nsample, length(gsvd.S), param_redu.DoF, d, md]);
        
        if md < redu_def.gsvd_conv_tol && gsvd.Nsample > redu_def.gsvd_Nmin
            stop_iter   = ii;
            break;
        end
        
        if gsvd.Nsample > redu_def.gsvd_Nmax
            stop_iter   = ii;
            break;
        end
    end
    
end

out.param_redu = param_redu;
out.param_dist = param_dist;
out.num_basis  = num_basis;
out.stop_iter  = stop_iter;

if redu_def.debug_flag 
    out.gsvd  = gsvd;
end

disp('done adaptation')

if redu_def.sub_chain
    out.sub_cov.mu  = zeros(param_redu.DoF,1);
    out.sub_cov.C   = diag((param_redu.S.^2 + 1).^(-1));
    out.sub_cov.num = stop_iter;
else
    samples         = out.v_samples(1:out.j,:)*param_redu.P;
    out.sub_cov.mu  = mean(samples)';
    out.sub_cov.C   = cov(samples);
    out.sub_cov.num = size(samples,1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mcmc_def, redu_def, out] = process_input(mcmc_def, redu_def)

redu_def                    = lis_input(redu_def);
if ~isfield(redu_def, 'sub_chain')
    redu_def.sub_chain    = false;
end
%number of sample intervals for global svd
if ~isfield(redu_def, 'gsvd_Ninter')
    redu_def.gsvd_Ninter    = 200;
end
% initial burning steps
if ~isfield(redu_def, 'dyna_N')
    redu_def.dyna_N         = 5;
end

n                           = redu_def.gsvd_Nmax*2+(redu_def.dyna_N+1)*redu_def.dyna_N;
mcmc_def.nstep              = (n/2)*redu_def.gsvd_Ninter;

user.high_d                 = true;
% default proposal
if ~isfield(mcmc_def, 'proposal')
    mcmc_def.proposal       = 'Prior';
end
% default acceptance rate
if strcmp(mcmc_def.proposal, 'MALA')
    user.rate               = 0.58;
else
    user.rate               = 0.23;
end    
% default jump size for the null space part
if ~isfield(mcmc_def, 'dt')
    mcmc_def.dt             = 2;
end
% using reference point
if ~isfield(mcmc_def, 'ref_type')
    mcmc_def.ref_type       = 'init';
end
% using curvature only
if ~isfield(mcmc_def, 'use_curvature')
    mcmc_def.use_curvature  = true;
end

mcmc_def.projection         = [];
mcmc_def.save_all           = true;
[mcmc_def, out]             = mcmc_input(mcmc_def, user);

end
