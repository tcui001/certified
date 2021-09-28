function stat = build_cov(mcmc_def, stat)
% Compute the covariacen
% Tiangang Cui, 01/Otc/2013

% coompute the mean and covariance
stat.M = stat.sum/stat.num;
stat.C = stat.cross/(stat.num-1) - stat.M(:)*stat.M(:)';
%
if mcmc_def.use_curvature
    stat.d = stat.acurv;
    stat.V = eye(length(stat.d));
else
    [V, D]          = eig(stat.C);
    [stat.d, ind]   = sort(abs(diag(D)), 'descend');
     stat.V         = V(:,ind);
    % d = abs(d + 0.05*stat.acurv);
end

end