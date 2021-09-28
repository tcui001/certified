function prior = set_prior_GP(C, mu)

prior.type      = 'Dist';
prior.cov.C     = C;
prior.cov.P     = inv(C);
prior.cov.RC    = chol(C);
prior.DoF       = size(C,1);
prior.cov.type  = 'GP';
prior.mean_u    = mu;


prior.func.type = 'log';
prior.func.log_thres    = zeros(prior.DoF,1);

end
