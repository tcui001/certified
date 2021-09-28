function redu_def = lis_input(redu_def)
%PROCESS_INPUT_LIS
%
% Process the input arguments for building the LIS
%
% Tiangang Cui, 17/Jan/2014

% options for Hessian construction
% define if using full Hessian, or the Gauss-Newton Hessian

% initial guess
if ~isfield(redu_def, 'init')
    disp('Error: initial guess needed');
end
redu_def.np             = length(redu_def.init);
% max samples for the gsvd
if ~isfield(redu_def, 'gsvd_Nmax')
    redu_def.gsvd_Nmax  = 500;
end
% min sample
if ~isfield(redu_def, 'gsvd_Nmin')
    redu_def.gsvd_Nmin  = 100;
end

if redu_def.gsvd_Nmax < redu_def.gsvd_Nmin
    redu_def.gsvd_Nmin  = redu_def.gsvd_Nmax;
end

% local truncation threshold
if ~isfield(redu_def, 'eigen_tol')
    redu_def.eigen_tol  = 1E-2;
end
if ~isfield(redu_def, 'gsvd_trunc_tol')
    redu_def.gsvd_trunc_tol = sqrt(redu_def.eigen_tol);
end
if ~isfield(redu_def, 'eigen_Nmax')
    redu_def.eigen_Nmax = 100;
end
if ~isfield(redu_def, 'gsvd_conv_tol')
    redu_def.gsvd_conv_tol  = 1E-3;
end
%
if ~isfield(redu_def, 'debug_flag')
    redu_def.debug_flag = false;
end
%
if ~isfield(redu_def, 'method')
    redu_def.method     = 'SVD';
end
redu_def.gsvd_tol = 1E-5; %sqrt(1E-1);
if ~isfield(redu_def, 'gsvd_Nmax_basis')
    redu_def.gsvd_Nmax_basis = 600;
end
if ~isfield(redu_def, 'gsvd_Niter_basis')
    redu_def.gsvd_Niter_basis = 100;
end

end