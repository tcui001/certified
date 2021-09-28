function [gsvd, redu] = redu_reinit(redu_def, gsvd, sample)
% Reinitialize the reduced subspace
% Tiangang Cui, 01/Oct/2013

% first compute the local subspace, then determine if use the iterative
% update or the full update

switch redu_def.method
    case {'Eig'}
        [V_locp,d_locp] = redu_def.PPGNH(sample, redu_def.eigen_tol, redu_def.eigen_Nmax);
         s_locp         = sqrt(d_locp);
    case {'SVD'}
        [~, s_locp, V_locp] = redu_def.SVD(sample);
         jnd            = s_locp >= sqrt(redu_def.eigen_tol);
         s_locp         = s_locp(jnd);
         V_locp         = V_locp(:,jnd);
end

% set redu_def.iterative_update_number to redu_def.gsvd_Nmax if full update
% is used all the time

iter_SVD_flag = true;

if iter_SVD_flag
    gD              = gsvd.S(:).^2;
    T               = gsvd.V'*V_locp;
    [Q, R]          = qr(V_locp - gsvd.V*T, 0);
    tmp1            = scale_cols([T; R], s_locp);
    tmp2            = tmp1*tmp1' + diag([gD*gsvd.Nsample; zeros(length(s_locp), 1)]);
    [Phi, D]        = eig(tmp2); % this line may cause numreical instability
    [S_glo,ind]     = sort( sqrt(diag(D)) / sqrt(gsvd.Nsample + 1), 'descend' );
    V_glo           = [gsvd.V, Q]*Phi(:, ind);
else % full update
      Vmat          = [scale_cols(gsvd.V, gsvd.S)*sqrt(gsvd.Nsample), scale_cols(V_locp, s_locp)];
    [V_glo, S_glo]  = svd(Vmat/sqrt(gsvd.Nsample+1),'econ'); % global SVD
     S_glo          = diag(S_glo); 

end

ind             = S_glo>=redu_def.gsvd_tol;
if redu_def.gsvd_Niter_basis <= gsvd.Nsample
    ind( (redu_def.gsvd_Nmax_basis+1):end ) = false;
end
gsvd.V          = V_glo(:,ind);
gsvd.S          = S_glo(ind);
gsvd.Nsample    = gsvd.Nsample + 1;

%{
switch redu_def.method
    case {'Eig'}
        [V_locp, d_locp] = redu_def.PPGNH(sample, redu_def.eigen_tol, redu_def.eigen_Nmax);
        Vmat    = [scale_cols(gsvd.V, gsvd.S)*sqrt(gsvd.Nsample), scale_cols(V_locp, sqrt(d_locp))];
    case {'SVD'}
        [~, s_locp, V_locp] = redu_def.SVD(sample);
        jnd     = s_locp >= sqrt(redu_def.eigen_tol);
        Vmat    = [scale_cols(gsvd.V, gsvd.S)*sqrt(gsvd.Nsample), scale_cols(V_locp(:,jnd), s_locp(jnd))];
end
%}

ind_trunc           = gsvd.S >= redu_def.gsvd_trunc_tol;
redu.DoF            = sum(ind_trunc);
redu.P              = gsvd.V(:,ind_trunc);
redu.S              = gsvd.S(ind_trunc);

end

