function [curr, stat, kernel] = dili_init(mcmc_def, param_redu, sigma)
% Initialize the initialize the definition of proposals duirng the LIS update
% Tiangang Cui, 17/Jan/2014
%

% initialize the covariance in the subspace
stat.acurv  = (param_redu.S.^2 + 1).^(-1);

if isfield(param_redu, 'cov')
    stat.num    = param_redu.cov.num;
    mu          = param_redu.cov.mu;
    C           = param_redu.cov.C;
    stat.sum    = mu*stat.num;
    stat.cross  = (C + mu*mu')*(stat.num - 1);
    
else
    stat.num    = param_redu.DoF;
    stat.sum    = (param_redu.P'*mcmc_def.init) *stat.num;
    stat.cross  = diag( stat.acurv )*stat.num + stat.sum*stat.sum'/stat.num;
end

stat        = build_cov(mcmc_def, stat);
kernel      = build_kernel(mcmc_def, param_redu, stat, sigma);
% initialize the state
curr.v      = mcmc_def.init;

switch mcmc_def.proposal
    case {'MALA'}
        [curr.mlpt, curr.mllkd, curr.grad] = mcmc_def.minus_log_post(curr.v);
        curr.grad_v_sub = param_redu.P'*curr.grad;
    otherwise
        [curr.mlpt, curr.mllkd] = mcmc_def.minus_log_post(curr.v);
end

curr.v_sub  = param_redu.P'*curr.v;
curr.v_null = curr.v - param_redu.P*curr.v_sub;

%norm(curr.v)
%norm(curr.v_null)

end
