function Jty = RD_matvec_Jty(FEM, HI, dy)
%HEAT_MATVEC_JTY
%
% Adjoint solve
%
% Tiangang Cui, 09/May/2014 

N       = size(dy, 2)/FEM.Ndatasets;
Jty     = zeros(FEM.mesh.Nel, FEM.NSteps+1);
lambda  = zeros(FEM.DoF,      FEM.Nsteps+1);

for i   = 1:N
    ind = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    dU  = FEM.C'*(dy(:,ind)*HI.T');
    
    
    for j = (soln.Nend+1):-1:1
        if  j   == (soln.Nend+1)
            g   = - dU(:,j);
        else
            switch soln.solver_type
                case{'MG', 'MD'}
                    g   = FEM.M*lambda(:,j) - dU(:,j);
                case {'MLG', 'MLD'}
                    g   = lambda(:,j) - dU(:,j);
            end
        end
        
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
                ds      = gmres(@(ds) ML_matvec_dFdu(FEM, s, ds),   g, RESTART, TOL, MAXIT, @(w) precond(HI.Us{j-1}', w), @(w) precond(HI.Ls{j-1}', w));
            case {'MLD', 'MD'}
                lambda(:,j)     = HI.Ls{j-1}'\(HI.Us{j-1}'\g);
        end
        
        lambda(FEM.p,j) = HI.R{j}\(HI.R{j}'\g(FEM.p));
    end
    for j = 2:(FEM.Nsteps+1)
        T1      = reshape(lambda(FEM.mesh.node_map,j),4,FEM.mesh.Nel);
        T2      = reshape(  HI.G(FEM.mesh.node_map,j),4,FEM.mesh.Nel);
        Jty(:,i)= Jty(:,i) + sum(T1.*(FEM.mesh.locstiff*T2))'*FEM.Tsteps(j-1);
    end
end

end

