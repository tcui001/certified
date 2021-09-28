function soln = RD_solve(FEM, s_init)
%HEAT_SOLVE  
%
% solves the transient heat equation by Cholesky and implicit Euler
%
% Tiangang Cui, 09/May/2014

RESTART             = 3;
TOL                 = 1E-6;
MAXIT               = 50;
NR_tol              = 1E-3;

soln.G              = zeros(FEM.DoF, FEM.Nmaxsteps+1);
soln.G(:,1)         = s_init;

dt                  = FEM.dt_init;
time                = 0;
soln.Nend           = 0;
soln.dts            = zeros(FEM.Nmaxsteps, 1);

%solver_type        = 'MLG';
if FEM.ML_flag
    if FEM.GMRES_flag
        soln.solver_type    = 'MLG';
    else
        soln.solver_type    = 'MLD';
    end
else
    if FEM.GMRES_flag
        soln.solver_type    = 'MG';
    else
        soln.solver_type    = 'MD';
    end
end

if FEM.fast_adj_flag
    soln.Ls         = cell{FEM.Nmaxsteps, 1};
    soln.Us         = cell{FEM.Nmaxsteps, 1};
end

for i = 1:FEM.Nmaxsteps
    % Newton iteration
    NR_iter         = 0;
    s               = soln.G(:,i);
    while 1
        switch soln.solver_type
            case{'MG'}
                F       = full_eval_F(FEM, soln.G(:,i), s);
                [L,U]   = full_calc_precond(FEM, s);
                ds      = gmres(@(ds) full_matvec_dFdu(FEM, s, ds), -F, RESTART, TOL, MAXIT, @(w) precond(L, w), @(w) precond(U, w));
            case {'MD'}
                [ds,L,U]= full_direct_solve(FEM, s0, s, dt);
            case {'MLG'}
                F       = ML_eval_F(FEM, soln.G(:,i), s);
                [L,U]   = ML_calc_precond(FEM, s);
                ds      = gmres(@(ds) ML_matvec_dFdu(FEM, s, ds),   -F, RESTART, TOL, MAXIT, @(w) precond(L, w), @(w) precond(U, w));
            case {'MLD'}
                [ds,L,U]= ML_direct_solve(FEM, s0, s, dt);
        end
        s       = s+ds;
        NR_iter = NR_iter + 1;
        if norm(ds)/FEM.DoF < NR_tol
            break;
        end
    end
    
    if FEM.fast_adj_flag
        soln.Ls{i}  = L;
        soln.Us{i}  = U;
    end
    
    soln.G(:,i+1)   = s; 
    time            = time + dt;
    soln.dts(i)     = dt;
    
    if time >= FEM.Tfinal
        soln.Nend   = i;
        break;
    end
    
    if NR_iter <= 2
        dt          = dt * FEM.dt_multiplier;
    end
end

% process observations 
obs_Tsteps          = linspace(0, FEM.Tfinal, FEM.Ndatasets);
soln.T              = sparse([],[],[], soln.Nend+1, FEM.Ndatasets, FEM.Ndatasets*2);
t_start             = 0;
for i = 1:soln.Nend
    t_end           = t_start + soln.dts(i);
    t_ind           = find(obs_Tsteps>t_start & obs_Tsteps<=t_end);
    temp1           = (t_end  -  obs_Tsteps(t_ind))/soln.dts(i); % weighting at t_start
    temp2           = (obs_Tsteps(t_ind) - t_start)/soln.dts(i); % weighting at t_end
    t_start         = t_end; % increment the time
    
    soln.T(i,  t_ind)   = temp1;
    soln.T(i+1,t_ind)   = temp2;    
end

soln.G              = soln.G(:,1:soln.Nend+1);
soln.d              = (FEM.C*soln.G)*soln.T;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w      = full_eval_F(FEM, s0, s, dt)
    w               = FEM.M*s - FEM.M*(s0) + dt*(FEM.K*s) - dt*FEM.a*(FEM.M*(s.*(1-s)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, U] = full_calc_precond(FEM, s, dt)
    % dFdu          = (FEM.M + FEM.K*dt) - FEM.a*FEM.M + 2*FEM.a*FEM.M*diag(s)
    S               = spdiags(s(:), 0, FEM.DoF, FEM.DoF);
    J               = FEM.M + dt*FEM.K - dt*FEM.a*FEM.M + 2*dt*FEM.a*FEM.M*S;
    setup.type      = 'crout';
    setup.milu      = 'row';
    setup.droptol   = 0.1;
    [L, U]          = ilu(J, setup);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w      = full_matvec_dFdu(FEM, s, ds, dt)
    % dFdu          = (FEM.M + FEM.K*dt) - FEM.a*FEM.M + 2*FEM.a*FEM.M*diag(s)
    w               = FEM.M*ds + dt*(FEM.K*ds) - dt*FEM.a*(FEM.M*ds) + 2*dt*FEM.a*(FEM.M*(s.*ds));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ds, L, U]  = full_direct_solve(FEM, s0, s, dt)
    F               = full_eval_F(FEM, s0, s, dt);
    S               = spdiags(s(:), 0, FEM.DoF, FEM.DoF);
    J               = FEM.M + dt*FEM.K - dt*FEM.a*FEM.M + 2*dt*FEM.a*FEM.M*S;
    [L, U]          = lu(J);
    ds              = -U\(L\F);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w      = ML_eval_F(FEM, s0, s, dt)
    w               = s - s0 + dt*(FEM.iMLK*s) - dt*FEM.a*(s.*(1-s));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, U] = ML_calc_precond(FEM, s, dt)
    % dFdu          = (FEM.M + FEM.K*dt) - FEM.a*FEM.M + 2*FEM.a*FEM.M*diag(s)
    S               = spdiags(s(:), 0, FEM.DoF, FEM.DoF);
    %I              = speye(FEM.DoF);
    J               = FEM.I + dt*FEM.iMLK - dt*FEM.a*FEM.I + 2*dt*FEM.a*S;
    setup.type      = 'crout';
    setup.milu      = 'row';
    setup.droptol   = 0.1;
    [L, U]          = ilu(J, setup);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w      = ML_matvec_dFdu(FEM, s, ds, dt)
    % dFdu          = (FEM.M + FEM.K*dt) - FEM.a*FEM.M + 2*FEM.a*FEM.M*diag(s)
    w               = ds + dt*(FEM.iMLK*ds) - dt*FEM.a*ds + 2*dt*FEM.a*(s.*ds);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ds, L, U]  = ML_direct_solve(FEM, s0, s, dt)
    F               = ML_eval_F(FEM, s0, s, dt);
    S               = spdiags(s(:), 0, FEM.DoF, FEM.DoF);
    J               = FEM.I + dt*FEM.iMLK - dt*FEM.a*FEM.I + 2*dt*FEM.a*S;
    [L, U]          = lu(J);
    ds              = -U\(L\F);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function invMw  = precond(M, w)
    invMw           = M\w;
end
