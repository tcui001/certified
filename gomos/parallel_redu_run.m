% parallel_redu_run
% rename input to P_max
% output is P_mllkd

ts = zeros(prior.DoF, nrp);
for j = 1:nrp
    ind = 1:nps(j);
    ts(:,j) = r - P_max(:,ind)*(P_max(:,ind)'*r);
end
    
tic
spmd
    Ml = getLocalPart(Md);
    nl = size(Ml, 2);
    mllkd   = zeros(nrp, nl);
    
    % projected samples
    rpl = P_max'*Ml;
    
    for j = 1:nrp
        ind = 1:nps(j);
        for i = 1:nl    
            c2p = pre_process(prior, P_max(:,ind)*rpl(ind,i) + ts(:,j), false);
            
            HI          = forward_solve(model, c2p.x);
            misfit      = (HI.d - obs.data)./obs.std;
            mllkd(j,i)  = 0.5*sum(misfit(:).^2); % minus log-likelihood
        end
    end
    par_mllkd = mllkd;
end
toc

