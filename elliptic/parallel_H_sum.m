% Md samples need to be distributed

tic
spmd
    Hl_g = zeros(prior.DoF);
    Hl_h = zeros(prior.DoF);
    Ml = getLocalPart(Md);
    
    for i = 1:size(Ml, 2)
        [~,~,~,gmllkd,HI] = minus_log_post(model, obs, prior, Ml(:,i));
        
        Hl_g = Hl_g + gmllkd*gmllkd';
        
        [V, D]      = myeig(HI);
        d           = diag(D);
        ind         = d>1E-6;
        
        Hl_h = Hl_h + V(:,ind)*scale_cols(V(:,ind),d(ind))';
    end
    par_Hg = Hl_g;
    par_Hh = Hl_h;
end
toc

tmp_Hg = zeros(prior.DoF);
tmp_Hh = zeros(prior.DoF);
tic;
for i = 1:nlabs
    tmp_Hg = tmp_Hg + par_Hg{i};
    tmp_Hh = tmp_Hh + par_Hh{i};
end
toc

