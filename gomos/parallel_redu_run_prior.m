% parallel_redu_run
% rename input to P_max
% output is P_mllkd

pr_mllkd        = zeros(nrp, out_dili.size);

prior_KLs = cell(nrp, 1);

ts = zeros(prior.DoF, nrp);
for j = 1:nrp
    prior_KLs{j} = basis_KL(prior, nps(j));
    ind = 1:nps(j);
    ts(:,j) = ur - prior_KLs{j}.P*(prior_KLs{j}.P'*ur) + prior.mean_u;
end
    
tic
spmd
    Ml = getLocalPart(Md);
    nl = size(Ml, 2);
    mllkd   = zeros(nrp, nl);
    uMl = matvec_prior_L(prior, Ml);
    
    % projected samples
    
    for j = 1:nrp
        rpl = prior_KLs{j}.P'*uMl;
        for i = 1:nl    
            u = prior_KLs{j}.P*rpl(:,i) + ts(:,j);
            x = u2x( prior, u );
            
            HI          = forward_solve(model, x);
            misfit      = (HI.d - obs.data)./obs.std;
            mllkd(j,i)  = 0.5*sum(misfit(:).^2); % minus log-likelihood
        end
    end
    par_mllkd = mllkd;
end
toc


for i = 1:nlabs
    ind = (1:seg) + (i-1)*seg;
    pr_mllkd(:,ind)    = par_mllkd{i};
end