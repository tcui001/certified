
n_trial         = 10;

pr_mllkds       = cell(n_trial,1);
new_mllkds      = cell(n_trial,1);
lis_mllkds      = cell(n_trial,1);
as_mllkds       = cell(n_trial,1);
prlis_mllkds    = cell(n_trial,1);
lap_mllkds      = cell(n_trial,1);
laplis_mllkds   = cell(n_trial,1);

load('mllkd_data_full.mat')
for kk = 1:n_trial
    fstr = ['mllkd_data_rand_' num2str(kk), '.mat'];
    load(fstr)
    pr_mllkds{kk}       = pr_mllkd;
    new_mllkds{kk}      = new_mllkd;
    lis_mllkds{kk}      = lis_mllkd;
    as_mllkds{kk}       = as_mllkd;
    prlis_mllkds{kk}    = prlis_mllkd;
    lap_mllkds{kk}      = lap_mllkd;
    laplis_mllkds{kk}   = laplis_mllkd;
end

nps = 10:5:50;
nrp = length(nps);

pr_dkl      = zeros(n_trial, nrp);
pr_dh       = zeros(n_trial, nrp);
new_dkl     = zeros(n_trial, nrp);
new_dh      = zeros(n_trial, nrp);
lis_dkl     = zeros(n_trial, nrp);
lis_dh      = zeros(n_trial, nrp);
as_dkl      = zeros(n_trial, nrp);
as_dh       = zeros(n_trial, nrp);
prlis_dkl   = zeros(n_trial, nrp);
prlis_dh    = zeros(n_trial, nrp);
lap_dkl     = zeros(n_trial, nrp);
lap_dh      = zeros(n_trial, nrp);
laplis_dkl  = zeros(n_trial, nrp);
laplis_dh   = zeros(n_trial, nrp);

min_pr_mllkds       = ones(nrp*out_dili.size, 1)*1E7;
min_new_mllkds      = ones(nrp*out_dili.size, 1)*1E7;
min_lis_mllkds      = ones(nrp*out_dili.size, 1)*1E7;
min_as_mllkds       = ones(nrp*out_dili.size, 1)*1E7;
min_prlis_mllkds    = ones(nrp*out_dili.size, 1)*1E7;
min_lap_mllkds      = ones(nrp*out_dili.size, 1)*1E7;
min_laplis_mllkds   = ones(nrp*out_dili.size, 1)*1E7;
for kk = 1:n_trial
    [pr_dkl(kk,:), pr_dh(kk,:)]         = calc_distances(full_mllkd, pr_mllkds{kk});
    [new_dkl(kk,:), new_dh(kk,:)]       = calc_distances(full_mllkd, new_mllkds{kk});
    [lis_dkl(kk,:), lis_dh(kk,:)]       = calc_distances(full_mllkd, lis_mllkds{kk});
    [as_dkl(kk,:), as_dh(kk,:)]         = calc_distances(full_mllkd, as_mllkds{kk});
    [prlis_dkl(kk,:), prlis_dh(kk,:)]   = calc_distances(full_mllkd, prlis_mllkds{kk});
    [lap_dkl(kk,:), lap_dh(kk,:)]       = calc_distances(full_mllkd, lap_mllkds{kk});
    [laplis_dkl(kk,:), laplis_dh(kk,:)] = calc_distances(full_mllkd, laplis_mllkds{kk});
    
    min_pr_mllkds       = min([pr_mllkds{kk}(:), min_pr_mllkds(:)], [], 2);
    min_new_mllkds      = min([new_mllkds{kk}(:), min_new_mllkds(:)], [], 2);
    min_lis_mllkds      = min([lis_mllkds{kk}(:), min_lis_mllkds(:)], [], 2);
    min_as_mllkds       = min([as_mllkds{kk}(:), min_as_mllkds(:)], [], 2);
    min_prlis_mllkds    = min([prlis_mllkds{kk}(:), min_prlis_mllkds(:)], [], 2);
    min_lap_mllkds      = min([lap_mllkds{kk}(:), min_lap_mllkds(:)], [], 2);
    min_laplis_mllkds   = min([laplis_mllkds{kk}(:), min_laplis_mllkds(:)], [], 2);
end

min_pr_mllkds       = reshape(min_pr_mllkds, nrp, out_dili.size);
min_new_mllkds      = reshape(min_new_mllkds, nrp, out_dili.size);
min_lis_mllkds      = reshape(min_lis_mllkds, nrp, out_dili.size);
min_as_mllkds       = reshape(min_as_mllkds, nrp, out_dili.size);
min_prlis_mllkds    = reshape(min_prlis_mllkds, nrp, out_dili.size);
min_lap_mllkds      = reshape(min_lap_mllkds, nrp, out_dili.size);
min_laplis_mllkds   = reshape(min_laplis_mllkds, nrp, out_dili.size);

% expected llkd 
pr_mllkde       = zeros(nrp, out_dili.size);
new_mllkde      = zeros(nrp, out_dili.size);
lis_mllkde      = zeros(nrp, out_dili.size);
as_mllkde       = zeros(nrp, out_dili.size);
prlis_mllkde    = zeros(nrp, out_dili.size);
lap_mllkde      = zeros(nrp, out_dili.size);
laplis_mllkde   = zeros(nrp, out_dili.size);
for kk = 1:n_trial
    pr_mllkde       = pr_mllkde + exp(min_pr_mllkds - pr_mllkds{kk});
    new_mllkde      = new_mllkde + exp(min_new_mllkds - new_mllkds{kk});
    lis_mllkde      = lis_mllkde + exp(min_lis_mllkds - lis_mllkds{kk});
    as_mllkde       = as_mllkde + exp(min_as_mllkds - as_mllkds{kk});
    prlis_mllkde    = prlis_mllkde + exp(min_prlis_mllkds - prlis_mllkds{kk});
    lap_mllkde      = lap_mllkde + exp(min_lap_mllkds - lap_mllkds{kk});
    laplis_mllkde   = laplis_mllkde + exp(min_laplis_mllkds - laplis_mllkds{kk});
end
pr_mllkde       = - log(pr_mllkde/n_trial) + min_pr_mllkds;
new_mllkde      = - log(new_mllkde/n_trial) + min_new_mllkds;
lis_mllkde      = - log(lis_mllkde/n_trial) + min_lis_mllkds;
as_mllkde       = - log(as_mllkde/n_trial) + min_as_mllkds;
prlis_mllkde    = - log(prlis_mllkde/n_trial) + min_prlis_mllkds;
lap_mllkde      = - log(lap_mllkde/n_trial) + min_lap_mllkds;
laplis_mllkde   = - log(laplis_mllkde/n_trial) + min_laplis_mllkds;


[pr_dkle, pr_dhe]           = calc_distances(full_mllkd, pr_mllkde);
[new_dkle, new_dhe]         = calc_distances(full_mllkd, new_mllkde);
[lis_dkle, lis_dhe]         = calc_distances(full_mllkd, lis_mllkde);
[as_dkle, as_dhe]           = calc_distances(full_mllkd, as_mllkde);
[prlis_dkle, prlis_dhe]     = calc_distances(full_mllkd, prlis_mllkde);
[lap_dkle, lap_dhe]         = calc_distances(full_mllkd, lap_mllkde);
[laplis_dkle, laplis_dhe]   = calc_distances(full_mllkd, laplis_mllkde);

% zero llkd

fstr = 'mllkd_data_fixed.mat';
load(fstr)
[pr_dklf, pr_dhf]           = calc_distances(full_mllkd, pr_mllkd);
[new_dklf, new_dhf]         = calc_distances(full_mllkd, new_mllkd);
[lis_dklf, lis_dhf]         = calc_distances(full_mllkd, lis_mllkd);
[as_dklf, as_dhf]           = calc_distances(full_mllkd, as_mllkd);
[prlis_dklf, prlis_dhf]     = calc_distances(full_mllkd, prlis_mllkd);
[lap_dklf, lap_dhf]         = calc_distances(full_mllkd, lap_mllkd);
[laplis_dklf, laplis_dhf]   = calc_distances(full_mllkd, laplis_mllkd);
    


% load data for computing the error bound
%load H_data2.mat

%nps = 10:5:60;

%dd = diag(new_D);
%new_bound = flipud(cumsum(dd(end:-1:1)));

pr_bounds = calc_bounds_pr(new_H, prior, nps);
new_bounds = calc_bounds(new_H, new_H, nps);
lis_bounds = calc_bounds(new_H, lis_H, nps);
as_bounds = calc_bounds(new_H, as_H, nps);
prlis_bounds = calc_bounds(new_H, prlis_H, nps);
lap_bounds = calc_bounds(new_H, lap_H, nps);
laplis_bounds = calc_bounds(new_H, laplis_H, nps);
