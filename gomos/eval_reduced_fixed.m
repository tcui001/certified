

nlabs = 32;
parpool('local', nlabs)

Md  = distributed(out_dili.v_samples');
par_mllkd   = Composite(nlabs);


% dimensionalities used in testing
nps = 10:5:50;
nrp = length(nps);

% output data
full_mllkd      = zeros(1,   out_dili.size);
pr_mllkd        = zeros(nrp, out_dili.size);
new_mllkd       = zeros(nrp, out_dili.size);
lis_mllkd       = zeros(nrp, out_dili.size);
as_mllkd        = zeros(nrp, out_dili.size);
prlis_mllkd     = zeros(nrp, out_dili.size);
lap_mllkd       = zeros(nrp, out_dili.size);
laplis_mllkd    = zeros(nrp, out_dili.size);

seg             = out_dili.size/nlabs;

build_redu_params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full llkd
%{
tic
spmd
    Ml = getLocalPart(Md);
    nl = size(Ml, 2);
    
    mllkd   = zeros(1, nl);
    
    for i = 1:nl
        c2p = pre_process(prior, Ml(:,i), false);
        HI          = forward_solve(model, c2p.x);
        misfit      = (HI.d - obs.data)./obs.std;
        mllkd(i)    = 0.5*sum(misfit(:).^2); % minus log-likelihood
    end
    
    par_mllkd = mllkd;
end
toc

for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    full_mllkd(ind)     = par_mllkd{i};
end

save mllkd_data_full.mat full_mllkd
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parallel_redu_run_prior_fixed;

% H matrix
P_max = new_P;
 parallel_redu_run_fixed;
for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    new_mllkd(:,ind)    = par_mllkd{i};
end

% H_lis matrix
P_max = lis_P;
 parallel_redu_run_fixed;
for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    lis_mllkd(:,ind)    = par_mllkd{i};
end

% H_as matrix
P_max = as_P;
 parallel_redu_run_fixed;
for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    as_mllkd(:,ind)     = par_mllkd{i};
end

% H_prlis matrix
P_max = prlis_P;
 parallel_redu_run_fixed;
for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    prlis_mllkd(:,ind)  = par_mllkd{i};
end

% H_lap matrix
P_max = lap_P;
 parallel_redu_run_fixed;
for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    lap_mllkd(:,ind)    = par_mllkd{i};
end

% H_laplis matrix
P_max = laplis_P;
 parallel_redu_run_fixed;
for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    laplis_mllkd(:,ind) = par_mllkd{i};
end

fstr = ['mllkd_data_fixed.mat'];
save(fstr, 'r', 'pr_mllkd', 'new_mllkd', 'lis_mllkd', 'as_mllkd',  'prlis_mllkd', 'lap_mllkd', 'laplis_mllkd' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%