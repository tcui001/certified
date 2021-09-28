function out = pcn_mcmc(mcmc_def)
%PCN_MCMC 
%
% Prior preconditioned Crank-Nicolson MCMC
% Tiangang Cui, 5/Mar/2013
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
% proposal:
%       Type of proposals, includes 'MALA', 'RW', default is 'RW'.
%
% rate:
%       Target accceptance rate default is 0.23.
%
% sigma:
%       To set the initial dt 
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
% save_all_flag:
%       A flag indicates save the entire MCMC trace of parameters, default
%       is false. Only turn on with sufficient amount of memory
%
% Pass-in function required:
%
% [minus_log_post, minus_log_likelihood, grad] = minus_log_post(v)
%
% In the examples we have, the following is used
% minus_log_post = @(v) minus_log_post_mcmc(mesh, obs, prior, param, FEM, v)
%
% input:
%       v, parameters equipped with identity prior covariance
%
% output:
%       minus_log_post, minus_log_likelihood, and grad, where grad is the
%       gradient of minus_log_post w.r.t. v
%
% This function is passed-in by the structure mcmc_def, e.g.,
% mcmc_def.minus_log_post = @(v) minus_log_post_mcmc(..., v);
% where ... are used provided data structures for setup the forward
% simulation, and evaluating posterior distribution, this function only
% pass-in one variable v.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[def, out]  = process_input(mcmc_def);
sigma       = def.sigma;
curr.v      = def.init;

% initialize the MCMC
switch def.proposal
    case {'RW'}
        [curr.mlpt, curr.mllkd] = def.minus_log_post(curr.v);
    case {'MALA'}
        [curr.mlpt, curr.mllkd, curr.grad] = def.minus_log_post(curr.v);
        curr.grad   = curr.grad - curr.v;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCMC
acc         = 0;
batch       = 0;
out.j       = 0; % Added by Gianluca
for i = 1:def.nstep
    % propose, and evaluate the acceptance rate
    [alpha, next]   = propose(def, sigma, curr);
    if  log(rand)   < alpha
        curr        = next;
        acc         = acc+1;
    end
    batch           = batch + 1;
    nbatch          = 50; % Gianluca added this
    if  batch       == nbatch % Where is nbatch defined?
        delta       = min(0.1,sqrt(nbatch/i)); % adjust jump size
        if (acc/nbatch) < rate
            sigma   = sigma - delta;
        else
            sigma   = sigma + delta;
        end
        batch               = 0;
        acc                 = 0;
        out.k               = out.k+1;
        out.sigma(out.k)    = sigma;
        out.acc(out.k,1)    = acc/nbatch;
    end
    
    % save    
    if def.save_all && mod(i, def.save_batch) == 0
        out.j                   = out.j + 1;
        out.v_samples(out.j,:)  = v_curr;
    end
    if def.out_size > 0
        out.p_samples(i,:)  = def.projection'*v_curr;
    end
    out.llkd(i)             = mllkd_curr;
    out.lpt (i)             = mlpt_curr;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose(def, sigma, curr)

dt  = exp(sigma);
a   = 2*sqrt(2*dt)/(2+dt); % compute beta
b   = (2-dt)/(2+dt); % sqrt( 1-beta^2 )
c   = -2*dt/(2+dt);  % sqrt( 1-beta^2 ) - 1

% propose a condidate, then evaluate the MH ratio
switch def.proposal
    case {'RW'}
        r           = randn(length(curr.v),1);
        next.v      = b*curr.v + a*r;
        % next.lpotential = def.lpotential(next.v);
        [next.mlpt, next.mllkd] = def.minus_log_post(next.v);
        alpha       = curr.mllkd - next.mllkd;
    case {'MALA'}
        r           = randn(length(curr.v),1);
        next.v      = b*curr.v + c*curr.grad + a*r;   
        [next.mlpt, next.mllkd, next.grad] = def.minus_log_post(next.v);
        next.grad   = next.grad - next.v;
        
        % reverse
        log_yx      = - next.mllkd - 0.25*dt*norm(next.grad)^2 ...
                      - 0.5*next.grad(:)'*(curr.v(:) - next.v(:)) ...
                      - 0.25*dt*next.grad(:)'*(curr.v(:) + next.v(:));
        
        % forward
        log_xy      = - curr.mllkd - 0.25*dt*norm(curr.grad)^2 ...
                      - 0.5*curr.grad(:)'*(next.v(:) - curr.v(:)) ...
                      - 0.25*dt*curr.grad(:)'*(next.v(:) + curr.v(:));
        
        alpha       = log_yx - log_xy;
        
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [def, out] = process_input(def)

% default proposal
if ~isfield(def, 'proposal')
    def.proposal    = 'RW';
end

user.rate           = 0.23;
%user.sigma          = log(1/sqrt(np)); % What is np?
user.sigma          = 1.25;
user.high_d         = false;
[def, out]          = mcmc_input(def, user);

end