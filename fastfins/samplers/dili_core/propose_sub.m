function [alpha, next] = propose_sub(mcmc_def, param_redu, kernel, curr, r_sub)
% subspace proposals for Hessian separated MCMC
% Tiangang Cui, 25/mar/2013

switch mcmc_def.proposal
    case {'MALA'}
        [alpha, next] = propose_sub_mala(mcmc_def, param_redu, kernel, curr, r_sub);
    case {'RW'}
        [alpha, next] = propose_sub_rw(mcmc_def, param_redu, kernel, curr, r_sub);
    case {'Prior'}
        [alpha, next] = propose_sub_prior(mcmc_def, param_redu, kernel, curr, r_sub);
    case {'Post'}
        [alpha, next] = propose_sub_post(mcmc_def, param_redu, kernel, curr, r_sub);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_sub_mala(mcmc_def, param_redu, kernel, curr, r_sub)

% copy the y_null
next.v_null = curr.v_null;
% drift
% drift_curr = -0.5*(kernel.C*curr.grad_v_sub);
drift_curr  = -(kernel.C*curr.grad_v_sub);
% next v in subspace
% next.v_sub  = (curr.v_sub - kernel.ref) + kernel.ref + drift_curr + kernel.L*r_sub;
next.v_sub  = curr.v_sub + drift_curr + kernel.L*r_sub;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_sub_rw(mcmc_def, param_redu, kernel, curr, r_sub)

% copy the y_null
next.v_null = curr.v_null;
% next v in subspace
% next.v_sub  = (curr.v_sub-kernel.ref) + kernel.L*r_sub + kernel.ref;
next.v_sub  = curr.v_sub + kernel.L*r_sub;
% next v
next.v      = param_redu.P*next.v_sub + next.v_null;
% eval density
[next.mlpt, next.mllkd] = mcmc_def.minus_log_post(next.v);
% acceptance
alpha       = curr.mlpt - next.mlpt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_sub_prior(mcmc_def, param_redu, kernel, curr, r_sub)

% copy the y_null
next.v_null = curr.v_null;
% next v in subspace
next.v_sub  = kernel.A*(curr.v_sub - kernel.ref) + kernel.B*r_sub + kernel.ref;
% next v
next.v      = param_redu.P*next.v_sub + next.v_null; 
% next.v = param_redu.P*next.v_sub;
% eval density
[next.mlpt, next.mllkd] = mcmc_def.minus_log_post(next.v);
alpha       = curr.mllkd - next.mllkd + sum((curr.v_sub - next.v_sub).*kernel.ref);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, next] = propose_sub_post(mcmc_def, param_redu, kernel, curr, r_sub)

% copy the y_null
next.v_null = curr.v_null;
% next v in subspace
next.v_sub  = kernel.a*(curr.v_sub - kernel.ref) + kernel.B*r_sub + kernel.ref;
% next v
next.v      = param_redu.P*next.v_sub + next.v_null; 
% next.v = param_redu.P*next.v_sub;
% eval density
[next.mlpt, next.mllkd] = mcmc_def.minus_log_post(next.v);
% acceptance prob
% additional term in the acceptance probability
r_next      = ((1 + kernel.a)/kernel.b)*(kernel.D*(curr.v_sub - next.v_sub)) + r_sub;
tmp         = 0.5*( sum(curr.v_sub.^2) - sum(next.v_sub.^2) + sum(r_sub.^2) - sum(r_next.^2) );
alpha       = tmp + curr.mllkd - next.mllkd;


r_next_t    = kernel.b*kernel.D*curr.v_sub - kernel.a*r_sub - ...
              (kernel.a*(kernel.a+1)/kernel.b)*kernel.D*kernel.ref;
          
disp(norm(r_next_t - r_next));

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
