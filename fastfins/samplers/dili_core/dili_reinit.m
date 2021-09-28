function [curr, stat, kernel] = dili_reinit(mcmc_def, curr, param_redu, sigma, samples, lag)
% reinitialize the definition of proposals duirng the LIS update
% Tiangang Cui, 25/Mar/2013

switch mcmc_def.proposal
    case {'MALA'}
        curr.grad_v_sub = param_redu.P'*curr.grad;
end

curr.v_sub  = param_redu.P'*curr.v;
curr.v_null = curr.v - param_redu.P*curr.v_sub;

% initialize the covariance in the subspace
stat.acurv  = (param_redu.S.^2 + 1).^(-1);

if ~isempty(samples)
    samples     = samples*param_redu.P;
    nsample     = size(samples,1);
    % initialize
    stat.num    = nsample*lag;
    mu          = mean(samples)';
    C           = cov(samples);
    stat.sum    = mu*stat.num;
    stat.cross  = (C + mu*mu')*(stat.num - 1);
else    
    stat.num    = param_redu.DoF;
    stat.sum    = (param_redu.P'*mcmc_def.init) *stat.num;
    stat.cross  = diag( stat.acurv )*stat.num + stat.sum*stat.sum'/stat.num;
end

stat        = build_cov(mcmc_def, stat);
kernel      = build_kernel(mcmc_def, param_redu, stat, sigma);

end
