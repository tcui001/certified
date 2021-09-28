function out = dili_mcmc(mcmc_def, param_redu)
%DILI_MCMC
%
% Dimension invariant, likelihood informed (DILI) MCMC sampler
% Tiangang Cui, 01/Oct/2013
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
% Tuning parameters are given by mcmc_def structure, it has the following
% inputs
%
% nstep:
%       Number of MCMC steps, default is 10,000.
%
% nabtch:
%       Number of batch size used for adaptive changing the mcmc jumpsize
%       and empirical estimation of the covariance, default is 50.
%
% using_gibbs:
%       A binary flag to switch between all-at-once update and Metropolis-
%       within-Gibbs update, default is ``true''.
%
% proposal:
%       Type of proposals, includes 'Prior', 'Post', 'MALA', 'RW', default
%       is 'Prior'.
%       'Prior':
%         From preconditioned Langevin, using an empirical approximation to
%         the posterior covariance as precondition operator. This scheme
%         keeps the prior dusttribution as the invariant distribution
%       'Post' :
%         This scheme keeps the Gaussian distribution defined by the low
%         rank approximation of the posterior covairance as the invariant
%         distribution
%       'MALA' :
%         MALA in the subspace
%
% rate:
%       Target accceptance rate in the LIS, or the target acceptance rate
%       for the non-gibbs update, default is 0.23, 0.58 is set as default
%       for the gibbs MALA
%
% sigma:
%       To set the initial dt in the LIS
%
% dt:
%       To define the discretization level of the Langevin SDE in the CS,
%       default is 2, corresponding to the independece proposal in the CS
%
% ref_type:
%       Reference point, default is init
%       'prior_mean'    uses the prior mean,
%       'init'          uses the MAP point,
%       'post_mean'     uses the empirical posterior mean
%
% use_curvature:
%       Turn off the adaptation of the empirical covariance in the LIS,
%       default is false
%
% Options for saving outputs:
%
% save_bacth:
%       The batch size for saving mcmc history, default is 1.
%
% projection:
%       A set of vectors for projecting the MCMC sample trace of the
%       parameters, default is [], for memory saving.
%
% save_all:
%       A flag indicates save the entire MCMC trace of parameters, default
%       is false. Only turn on with sufficient amount of memory
%
% Pass-in function required:
%
% minus_log_post = @(y) minus_log_post_mcmc(mesh, obs, prior, param, FEM, y)
%
% [minus_log_post, minus_log_likelihood, grad] = minus_log_post(y)
%
% input:
%       y, parameters equipped with identity prior covariance
%
% output:
%       minus_log_post, minus_log_likelihood, and grad, where grad is the
%       gradient of minus_log_post w.r.t. y
%
% This function is passed-in by the strcuture mcmc_def, e.g.,
% mcmc_def.minus_log_post = @(y) minus_log_post_mcmc(..., y);
% where ... are used provided data structures for setup the forward
% simulation, and evaluating posterior distribution, this function only
% pass-in one variable y.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[def, out]  = process_input(mcmc_def); % parse inputs
sigma       = def.sigma;

[curr,stat,kernel]  = dili_init(def, param_redu, sigma); % initialization

% 
out.grad2   = zeros(size(param_redu.P, 1));

% MCMC
acc_sub     = 0;
acc_null    = 0;
batch       = 0;
%debug_flag
for i = 1:def.nstep
    % random number
    r                   = randn(def.np,1);
    
    if def.using_gibbs % gibbs update
        % project random number
        r_sub           = param_redu.P'*r;
        r_null          = r - param_redu.P*r_sub;
        
        % propose in the subspace, and evaluate the acceptance rate
        [alpha, next]   = propose_sub(def, param_redu, kernel, curr, r_sub);
        if  log(rand)   < alpha
            acc_sub     = acc_sub+1;
            curr        = next; % update
        end
        stat.cross      = stat.cross + curr.v_sub(:)*curr.v_sub(:)';
        stat.sum        = stat.sum + curr.v_sub;
        stat.num        = stat.num + 1;

        % propose in the null, and evaluate the acceptance rate
        [alpha, next]   = propose_null(def, param_redu, kernel, curr, r_null);
        if  log(rand)   < alpha
            acc_null    = acc_null+1;
            curr        = next; % update
        end
    else % all in one update
        % propose in the null, and evaluate the acceptance rate
        [alpha, next]   = propose(def, param_redu, kernel, curr, r);
        if  log(rand)   < alpha
            acc_sub     = acc_sub+1;
            curr        = next; % update
        end
        stat.cross      = stat.cross + curr.v_sub(:)*curr.v_sub(:)';
        stat.sum        = stat.sum + curr.v_sub;
        stat.num        = stat.num + 1;
    end % end of update
    
    out.grad2           = out.grad2 + curr.grad*curr.grad';
    
    %{
    % save
    if  def.save_batch  == 1
        out.j                       = out.j+1;
        if def.out_size > 0
            out.p_samples(out.j,:)  = def.projection'*curr.v;
        end
        if def.save_all_flag
            out.v_samples(out.j,:)  = curr.v;
        end
        out.llkd(out.j)             = curr.mllkd;
        out.lpt (out.j)             = curr.mlpt;
    elseif mod(i, def.save_batch)   == 0
        out.j                       = out.j+1;
        if def.out_size > 0
            out.p_samples(out.j,:)  = def.projection'*curr.v;
        end
        if def.save_all_flag
            out.v_samples(out.j,:)  = curr.v;
        end
        out.llkd(out.j)             = curr.mllkd;
        out.lpt (out.j)             = curr.mlpt;
    end
    %}
    
    % save    
    if def.save_all && mod(i, def.save_batch) == 0
        out.j                   = out.j + 1;
        out.v_samples(out.j,:)  = curr.v;
    end
    if def.out_size > 0
        out.p_samples(i,:)  = def.projection'*curr.v;
    end
    out.llkd(i)             = curr.mllkd;
    out.lpt (i)             = curr.mlpt;
    
    % tune the transition kernel
    batch       = batch + 1;
    if  batch   == def.nbatch
        delta   = min(0.1,sqrt(def.nbatch/i));
        if (acc_sub/def.nbatch) < def.rate
            sigma   = sigma - delta;
        else
            sigma   = sigma + delta;
        end
        
        % save
        out.k               = out.k+1;
        out.sigma(out.k,:)  = sigma;
        out.acc(out.k,:)    = [acc_sub, acc_null]/def.nbatch;
        
        % reset coounter
        batch       = 0;
        acc_sub     = 0;
        acc_null    = 0;
        % build new kernel, only when the adaptation is continue
        if i <= def.stop_adapt_n
            stat    = build_cov(def, stat);
        end
        kernel      = build_kernel(def, param_redu, stat, sigma);
    end
end

out.stat = stat;
out.grad2 = out.grad2/def.nstep;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [def, out] = process_input(def)

user.high_d         = true;
% default update methods, using_gibbs = true
if ~isfield(def, 'using_gibbs')
    def.using_gibbs = true;
end
% default proposal
if ~isfield(def, 'proposal')
    def.proposal    = 'Prior';
end
% default acceptance rate
if def.using_gibbs && strcmp(def.proposal, 'MALA')
    user.rate       = 0.58;
else
    user.rate       = 0.23;
end    
% default jump size for the null space part
if ~isfield(def, 'dt')
    def.dt          = 2;
end
% flag that indicate if adaptation on the covariance is turned off after
% certain number of iterations
if ~isfield(def, 'stop_adapt')
    def.stop_adapt  = true;
end
if def.stop_adapt
    def.stop_adapt_n    = def.nstep/4;
else
    def.stop_adapt_n    = def.nstep;
end

% using reference point
if ~isfield(def, 'ref_type')
    def.ref_type    = 'init';
end
% using curvature only
if ~isfield(def, 'use_curvature')
    def.use_curvature   = false;
end
[def, out]          = mcmc_input(def, user);

end
