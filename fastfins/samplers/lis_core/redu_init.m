function [gsvd, redu] = redu_init(redu_def, sample)
% Initialize the reduced subspace
% Tiangang Cui, 01/Oct/2013

% initialize

switch redu_def.method
    case {'Eig'}
        [V_loc, d_loc] = redu_def.PPGNH(sample, redu_def.eigen_tol, redu_def.eigen_Nmax);
        gsvd.V  = V_loc;
        gsvd.S  = d_loc.^(0.5);
    case {'SVD'}
        [~, s_loc, V_loc] = redu_def.SVD(sample);
        jnd     = s_locp >= sqrt(redu_def.eigen_tol);
        gsvd.V  = V_loc(:,jnd);
        gsvd.S  = s_loc(jnd);
end

gsvd.DoF        = length(gsvd.S);
gsvd.Nsample   = 1;
redu.DoF        = length(gsvd.S);
redu.P          = gsvd.V;
redu.S          = gsvd.S;

end