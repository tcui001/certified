function kernel = build_kernel(mcmc_def, param_redu, stat, sigma)
% Compute the operators for proposals duirng the LIS update
% Tiangang Cui, 25/Mar/2013
%

switch mcmc_def.proposal
    case {'MALA', 'RW'}
        kernel = build_kernel_sub  (stat.V, stat.d, sigma);
    case {'Prior'}
        kernel = build_kernel_prior(stat.V, stat.d, sigma);
    case {'Post'}
        kernel = build_kernel_post (stat.V, stat.d, sigma);
end

switch mcmc_def.ref_type
    case {'prior_mean'}
        kernel.ref = zeros(size(param_redu.P, 2), 1);
    case {'init'}
        kernel.ref = param_redu.P'*mcmc_def.init;   
    case {'post_mean'}
        kernel.ref = stat.M;
end

%{
switch mcmc_def.ref_type
    case {'prior_mean'}
        kernel.ref = zeros(size(mcmc_def.init));
        % kernel.ref_sub = zeros(size(stat.d));
    case {'init'}
        % kernel.ref = mcmc_def.init;
        kernel.ref = param_redu.P*(param_redu.P'*mcmc_def.init);   
    case {'post_mean'}
        %kernel.ref = mcmc_def.init;
        kernel.ref = param_redu.P*(stat.M+param_redu.P'*mcmc_def.init);  
        % kernel.ref = param_redu.P*M + ( mcmc_def.init - param_redu.P*(param_redu.P'*mcmc_def.init) );
        %kernel.ref_sub = stat.M;
end
%}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel = build_kernel_sub(V, d, sigma)

kernel.dt = exp(sigma)/sqrt(sum(d));
kernel.L  = scale_cols(V, d.^(0.5))*V'*sqrt(2*kernel.dt);
kernel.C  = scale_cols(V, d)*V'*kernel.dt;

%{
% for debug
scale = sqrt(2*kernel.dt);
kernel.dL = scale_cols(V, d.^(0.5))*V'*scale;
kernel.dC = scale_cols(V, d)*V'*scale^2;
%}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel = build_kernel_prior(V, d, sigma)
% build the kernel that have prior distribution as invariant, using
% preconditioned Langevin technique
% fix null space scaling, scale the kernel of the preconditioner P

kernel.dt = exp(sigma)/sqrt(sum(d));
tmp       = ( (0.5*kernel.dt)*d + 1 ).^(-1); % eigen values of the B operator
DB        = sqrt(2*kernel.dt)*sqrt(d).*tmp;
DA        = (1 - 0.5*kernel.dt*d).*tmp;
kernel.B  = scale_cols(V, DB)*V';
kernel.A  = scale_cols(V, DA)*V';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel = build_kernel_post(V, d, sigma)
% build the kernel that have Gaussian approximation of the posterior 
% distribution as invariant, using preconditioned Langevin technique
% fix null space scaling, scale the kernel of the preconditioner P

kernel.dt = exp(sigma)/sqrt(sum(d));
tmp       = 1/(1 + 0.5*kernel.dt); 
%DB        = sqrt(2*kernel.dt)*tmp*sqrt(d);
%kernel.B  = scale_cols(V, DB)*V';
kernel.a  = (1 - 0.5*kernel.dt)*tmp;
kernel.b  = sqrt(2*kernel.dt)*tmp;
kernel.D  = scale_cols(V, d.^(-0.5))*V'; % inv sqrt of the covariance
kernel.B  = scale_cols(V, kernel.b*sqrt(d))*V';

end
