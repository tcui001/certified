function redu_def = process_input_lis(redu_def)
%PROCESS_INPUT_LIS
%
% Process the input arguments for building the LIS
%
% Tiangang Cui, 17/Jan/2014

% options for Hessian construction
% define if using full Hessian, or the Gauss-Newton Hessian

redu_def.np = length(redu_def.init);
% max samples for the gsvd
if ~isfield(redu_def, 'gsvd_Nmax')
    redu_def.gsvd_Nmax = 500;
end
% min sample
if ~isfield(redu_def, 'gsvd_Nmin')
    redu_def.gsvd_Nmin = 100;
end
% local truncation threshold
if ~isfield(redu_def, 'eigen_tol')
    redu_def.eigen_tol = 1E-2;
end
redu_def.gsvd_trunc_tol = sqrt(redu_def.eigen_tol);
if ~isfield(redu_def, 'eigen_Nmax')
    redu_def.eigen_Nmax = 100;
end
if ~isfield(redu_def, 'gsvd_conv_tol')
    redu_def.gsvd_conv_tol = 1E-3;
end
%
if ~isfield(redu_def, 'debug_flag')
    redu_def.debug_flag = false;
end
%
if ~isfield(redu_def, 'method')
    redu_def.method = 'SVD';
end
redu_def.gsvd_tol = 1E-5; %sqrt(1E-1);
if ~isfield(redu_def, 'gsvd_Nmax_basis')
    redu_def.gsvd_Nmax_basis = 800;
end
if ~isfield(redu_def, 'gsvd_Niter_basis')
    redu_def.gsvd_Niter_basis = 20;
end

% algorithm specific options
if strcmp(redu_def.algo, 'RTO')% RTO options
    % using random walk, only used by RTO
    if ~isfield(redu_def, 'type')
        redu_def.type = 'MCMC';
    end
    redu_def.gsvd_Ninter = 1;
else % DILI options
    %number of sample intervals for global svd
    if ~isfield(redu_def, 'gsvd_Ninter')
        redu_def.gsvd_Ninter = 200;
    end
    % initial burning steps
    if ~isfield(redu_def, 'dyna_N')
        redu_def.dyna_N = 1;
    end
    redu_def.nstep = (redu_def.gsvd_Nmax+(redu_def.dyna_N+1)*redu_def.dyna_N/2)*redu_def.gsvd_Ninter;
    % using all the previous samples to readapt the
    if ~isfield(redu_def, 're_adapt')
        redu_def.re_adapt = true;
    end
end
%

end