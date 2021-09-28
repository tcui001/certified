function [alpha, next] = propose_null(mcmc_def, param_redu, kernel, curr, r_null)
% null space proposal for the Hessian separated MCMC
% Tiangang Cui, 25/Mar/2013
%
beta   = 2*sqrt(2*mcmc_def.dt)/(2 + mcmc_def.dt);
beta_p = (2 - mcmc_def.dt)/(2 + mcmc_def.dt);
% copy the y_sub
next.v_sub  = curr.v_sub;
% next y in null space, without ref point
next.v_null = beta_p*curr.v_null + beta*r_null;
%next.v_null = beta_p*(curr.v_null-kernel.ref_null) + beta*r_null + kernel.ref_null;
% next y
next.v      = param_redu.P*next.v_sub + next.v_null;
% eval density
switch mcmc_def.proposal
    case {'MALA'}
        [next.mlpt, next.mllkd, next.grad] = mcmc_def.minus_log_post(next.v);
        % project the gradient
        next.grad_v_sub = param_redu.P'*next.grad;
    otherwise
        [next.mlpt, next.mllkd] = mcmc_def.minus_log_post(next.v);
end
% acceptance
alpha = curr.mllkd - next.mllkd;

end
