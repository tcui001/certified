function [alpha, next] = propose(mcmc_def, param_redu, kernel, curr, r)
% subspace proposals for Hessian preconditioned MCMC
% Tiangang Cui, 08/August/2013

switch mcmc_def.proposal
    case {'MALA'}
        [alpha, next] = propose_mala (mcmc_def, param_redu, kernel, curr, r);
    case {'Prior'}
        [alpha, next] = propose_prior(mcmc_def, param_redu, kernel, curr, r);
    case {'Post'}
        [alpha, next] = propose_post (mcmc_def, param_redu, kernel, curr, r);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_prior(mcmc_def, param_redu, kernel, curr, r)

beta   = 2*sqrt(2*mcmc_def.dt)/(2 + mcmc_def.dt);
beta_p = (2 - mcmc_def.dt)/(2 + mcmc_def.dt);

% project random number
r_sub  = param_redu.P'*r;
r_null = r - param_redu.P*r_sub;

% next v in subspace
next.v_sub  = kernel.A*(curr.v_sub-kernel.ref) + kernel.ref + kernel.B*r_sub;
% next v in null space, without ref point
next.v_null = beta_p*curr.v_null + beta*r_null;
% next v
next.v      = param_redu.P*next.v_sub + next.v_null;

[next.mlpt, next.mllkd] = mcmc_def.minus_log_post(next.v);
% acceptance
alpha = curr.mllkd - next.mllkd + sum((curr.v_sub - next.v_sub).*kernel.ref);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_post(mcmc_def, param_redu, kernel, curr, r)

beta   = 2*sqrt(2*mcmc_def.dt)/(2 + mcmc_def.dt);
beta_p = (2 - mcmc_def.dt)/(2 + mcmc_def.dt);

% project random number
r_sub  = param_redu.P'*r;
r_null = r - param_redu.P*r_sub;

% next v in subspace
next.v_sub  = kernel.a*(curr.v_sub-kernel.ref) + kernel.ref + kernel.B*r_sub;
% next v in null space, without ref point
next.v_null = beta_p*curr.v_null + beta*r_null;
% next v
next.v      = param_redu.P*next.v_sub + next.v_null;

[next.mlpt, next.mllkd] = mcmc_def.minus_log_post(next.v);
r_next = (1+kernel.a)/kernel.b*(kernel.D*(curr.v_sub-next.v_sub)) + r_sub;
tmp    = 0.5*( sum(curr.v_sub.^2) - sum(next.v_sub.^2) + sum(r_sub.^2) - sum(r_next.^2) );
% acceptance prob
alpha  = tmp + curr.mllkd - next.mllkd; 

% debug
t1  = eye(size(kernel.B))*kernel.a + kernel.B - eye(size(kernel.B));
t2  = eye(size(kernel.B)) - eye(size(kernel.B))*kernel.a;
t3  = (kernel.B)^(-2);
w   = curr.v_sub;
wp  = next.v_sub;
wwp = - curr.mllkd - 0.5*sum((t3*w) .*(t1*w))  - 0.5*sum((t3*kernel.ref).*(t2*w));
wpw = - next.mllkd - 0.5*sum((t3*wp).*(t1*wp)) - 0.5*sum((t3*kernel.ref).*(t2*wp));

disp(alpha - (wpw - wwp));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_mala(mcmc_def, param_redu, kernel, curr, r)

beta   = 2*sqrt(2*mcmc_def.dt)/(2 + mcmc_def.dt);
beta_p = (2 - mcmc_def.dt)/(2 + mcmc_def.dt);

% project random number
r_sub  = param_redu.P'*r;
r_null = r - param_redu.P*r_sub;

% drift
drift_curr  = -(kernel.C*curr.grad_v_sub);
% next v in subspace
next.v_sub  = curr.v_sub + drift_curr + kernel.L*r_sub;
% next v in null space, without ref point
next.v_null = beta_p*curr.v_null + beta*r_null;
% next v
next.v      = param_redu.P*next.v_sub + next.v_null;
% eval density
[next.mlpt, next.mllkd, next.grad] = mcmc_def.minus_log_post(next.v);
% project the gradient
next.grad_v_sub = param_redu.P'*next.grad;
% backward random number
r_next = 0.5*kernel.L*(curr.grad_v_sub + next.grad_v_sub) - r_sub;
% additional term in the acceptance probability
tmp    = 0.5*( sum(curr.v_sub.^2) - sum(next.v_sub.^2) + sum(r_sub.^2) - sum(r_next.^2) );
% acceptance prob
alpha  = tmp + curr.mllkd - next.mllkd;

end
