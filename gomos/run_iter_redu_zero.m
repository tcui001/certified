% run iterative reduction using importance sampling

% change the maximum to 40

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iteration zero, prior or laplace

out_subs    = cell(1, n_iter);
Hs          = cell(1, n_iter+1);
prior_rs    = cell(1, n_iter);
mllkds      = cell(1, n_iter);
kk = 1;

tic;
gs  = zeros(prior.DoF, batch_size);
M   = randn(prior.DoF, batch_size);
for i = 1:batch_size
    [~,~,~,gmllkd,~] = minus_log_post(model, obs, prior, M(:,i));
    gs(:,i) = gmllkd;
end
gs = gs / sqrt(batch_size);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for i = 1:n_iter
    
    % build reduced prior
    [V, S, ~] = svd(gs);
    s = diag(S);
    j = find( cumsum(s.^2, 'reverse') < tol, 1);
    
    jm  = min(j, max_size);
    prior_r = basis_LIS(prior, V(:, 1:jm));
    
    disp([jm, j, max_size])
    
    % save
    Hs{kk,i} = gs;
    prior_rs{kk,i} = prior_r;
    
    % init MCMC
    vrmap               = get_map_matlab(model, obs, prior_r, randn(prior_r.DoF, 1));
    [~, s, V, ~]        = svd_explicit_WJ(model, obs, prior_r, vrmap, 1E-2, 0);
    
    sub_def.proposal    = 'MALA';
    sub_def.init        = vrmap;
    sub_def.sigma       = -1;
    sub_def.nstep       = nstep;
    
    sub_def.cov_def     = 'Restart';
    sub_def.mean        = sub_def.init;
    sub_def.sample_num  = prior_r.DoF*10;
    sub_def.cov         = inv(V*(diag(s.^2)*V') + eye(prior_r.DoF));
    sub_def.minus_log_post  = @(v) minus_log_post(model, obs, prior_r, v);
    
    % run MCMC
    tic;
    out_subs{kk,i}      = amcmc(sub_def);
    toc
    
    tic;
    ind         = randperm(size(out_subs{kk, i}.v_samples, 1), batch_size);
    sub_s       = out_subs{kk, i}.v_samples(ind,:)';
    rs          = randn(prior.DoF, batch_size);
    null_s      = rs - prior_r.P*(prior_r.P'*rs);
    appro_s     = prior_r.P*sub_s + null_s;
    f_mllkds    = zeros(batch_size, 1);
    for j = 1:batch_size
        [~,mllkd,~,gmllkd,~] = minus_log_post(model, obs, prior, appro_s(:,j));
        gs(:, j)    = gmllkd;
        f_mllkds(j) = mllkd;
    end
    log_ws      = out_subs{kk,i}.llkd(ind) - f_mllkds;
    mws         = max(log_ws);
    sws         = exp(mws) * sum( exp(log_ws - mws) );
    if std(exp(log_ws)) > 1E2 || ip_off
        disp('IP off')
        gs = gs / sqrt( batch_size );
    else
        gs = (gs.*repmat(exp(0.5*log_ws(:)'), prior.DoF, 1) ) / sqrt( sws );
    end
    toc
    
    
    % reduction error
    tic
    mllkd = zeros(out_dili.size, 1);
    rpl = prior_r.P'*(out_dili.v_samples');
    for k = 1:out_dili.size
        c2p = pre_process(prior_r, rpl(:,k), false);
        
        HI          = forward_solve(model, c2p.x);
        misfit      = (HI.d - obs.data)./obs.std;
        mllkd(k)  = 0.5*sum(misfit(:).^2); % minus log-likelihood
    end
    toc
    mllkds{kk, i} = mllkd;
    
end
Hs{kk,n_iter+1} = gs;


fstr = 'run_iter_zero';
if ip_off
    fstr = [fstr, '.mat'];
else
    fstr = [fstr, '_ip', '.mat'];
end
save(fstr, 'Hs', 'prior_rs', 'out_subs', 'mllkds', '-v7.3');