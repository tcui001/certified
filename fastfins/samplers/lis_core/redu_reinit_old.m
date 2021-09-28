function [gsvd, redu] = redu_reinit_old(redu_def, gsvd, sample)
% Reinitialize the reduced subspace
% Tiangang Cui, 01/Oct/2013

switch redu_def.method
    case {'Eig'}
        [V_locp, d_locp] = redu_def.PPGNH(sample, redu_def.eigen_tol, redu_def.eigen_Nmax);
        Vmat    = [scale_cols(gsvd.V, gsvd.S)*sqrt(gsvd.Nsample), scale_cols(V_locp, sqrt(d_locp))];
    case {'SVD'}
        [~, s_locp, V_locp] = redu_def.SVD(sample);
        jnd     = s_locp >= sqrt(redu_def.eigen_tol);
        Vmat    = [scale_cols(gsvd.V, gsvd.S)*sqrt(gsvd.Nsample), scale_cols(V_locp(:,jnd), s_locp(jnd))];
end

[V_glo, S_glo]  = svd(Vmat/sqrt(gsvd.Nsample+1),'econ'); % global SVD
 d_glo          = diag(S_glo); 
 ind            = (d_glo>=redu_def.gsvd_tol);
gsvd.V          = V_glo(:,ind);
gsvd.S          = d_glo(ind);
gsvd.Nsample   = gsvd.Nsample + 1;

ind_trunc       = gsvd.S >= redu_def.gsvd_trunc_tol;
redu.DoF        = sum(ind_trunc);
redu.P          = gsvd.V(:,ind_trunc);
redu.S          = gsvd.S(ind_trunc);

end

