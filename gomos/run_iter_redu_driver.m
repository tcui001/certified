
[model, obs, prior] = setup_gomos_GP();

n_trial     = 0;
n_iter      = 5;
max_size    = 40;
tol         = 1E-2;
batch_size  = 500;
nstep       = 2E4;

ip_off      = false;
run_iter_redu_zero;

max_size    = 30;
ip_off      = true;
run_iter_redu_zero;

plot_iter_zero

%{
% run 1
init_t      = 'prior';
ip_off      = true;
run_iter_redu;

% run 2
init_t      = 'prior';
ip_off      = false;
run_iter_redu;
%}

% run 3
%init_t      = 'laplace';
%ip_off      = true;
%run_iter_redu;

% run 4
%init_t      = 'laplace';
%ip_off      = false;
%run_iter_redu;

%{
% run 5
init_t      = 'prior';
run_iter_redu_ps;

% run 6
init_t      = 'laplace';
run_iter_redu_ps;
%}