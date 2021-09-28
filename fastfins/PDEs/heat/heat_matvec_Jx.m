function Jx = heat_matvec_Jx(FEM, HI, dx)
%HEAT_MATVEC_JX
%
% Adjoint solve
%
% Tiangang Cui, 09/May/2014 

TOL         = 1E-6;
MAXIT       = 50;

N   = size(dx, 2);
du  = zeros(FEM.DoF,      FEM.Nsteps+1);
Jx  = zeros(FEM.Nsensors, FEM.Ndatasets*N);

for i = 1:N
    dK      = sparse(1:size(dx,1),1:size(dx,1),dx(:,i));
    dN      = FEM.W1*dK*FEM.W1' + FEM.W2*dK*FEM.W2' + FEM.W3*dK*FEM.W3';
    
    du(:,1) = zeros(FEM.DoF,1);
    for j   = 1:FEM.Nsteps
        %   compute for j+1
        g   = FEM.M*du(:,j) - FEM.Tsteps(j)*dN*HI.G(:,j+1);
        
        if FEM.GMRES_flag
            B   = FEM.M + FEM.Tsteps(j)*HI.N;
            [du(:,j+1),~] = minres(B, g, TOL, MAXIT, HI.L{j+1}, HI.L{j+1}');
        else
            du(FEM.p,j+1) = HI.L{j+1}'\(HI.L{j+1}\g(FEM.p));
        end
        
    end
    
    ind     = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    Jx(:,ind)   = (FEM.C*du)*FEM.obs_int;
end

end

%{

TOL         = 1E-6;
MAXIT       = 50;

N           = size(dx, 2);
ds          = zeros(FEM.DoF*2, HI.Nend+1);
dRes        = zeros(FEM.DoF*2, HI.Nend+1);
Jx          = zeros(FEM.Nsensors, FEM.Ndatasets*N);


% dResdalpha
tmp_u       = HI.G(FEM.ind_u,:);
dResdalpha  = scale_cols(tmp_u.*(1-tmp_u), HI.dts(1:HI.Nend+1));
    
for i   = 1:N
    % multiply with purtubation
    dRes(:,:)               = 0;
    dRes(FEM.ind_u,2:end)   = scale_rows(dResdalpha(:,2:end), dx(:,i));
    
    ds(:,:)                 = 0;
    for j = 1:HI.Nend
        dt                  = HI.dts(j+1);
        mudt                = FEM.mu*dt;
        epdt                = FEM.epsilon*dt;
        adt                 = FEM.a*dt;
        bdt                 = FEM.b*dt;

        A                   = kron(diag([mudt, epdt]), FEM.iMLK) + kron([1, dt; -adt, 1+bdt], FEM.I);
        J                   = direct_dResds(FEM, HI.alpha, A, dt, HI.G(:,j+1));
        
        g                   = J'*( FEM.vecML.*( ds(:,i) + dRes(:,j+1) ) );
        
        if HI.pre_flag
            if FEM.GMRES_flag
                H           = J'*FEM.ML2*J; 
                ds(:,j+1)   = minres(H, g, TOL, MAXIT, HI.Ls{j+1}, HI.Ls{j+1}');
            else
                ds(:,j+1)   = HI.Ls{j+1}'\(HI.Ls{j+1}\g);
            end
        else
            H               = J'*FEM.ML2*J; 
            if FEM.GMRES_flag
                L           = make_precond_sym(H);
                ds(:,j+1)   = minres(H, g, TOL, MAXIT, L, L');
            else
                ds(:,j+1)   = H\g;
            end
        end
    end
    
    ind                     = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    Jx(:,ind)               = - (FEM.C*ds(FEM.ind_u,:))*HI.T;
end

%}


%{

N   = size(dx, 2);
du  = zeros(FEM.DoF,      FEM.Nsteps+1);
Jx  = zeros(FEM.Nsensors, FEM.Ndatasets*N);

for i = 1:N
    dK      = sparse(1:size(dx,1),1:size(dx,1),dx(:,i));
    dN      = FEM.W1*dK*FEM.W1' + FEM.W2*dK*FEM.W2' + FEM.W3*dK*FEM.W3';
    
    du(:,1) = zeros(FEM.DoF,1);
    for j   = 1:FEM.Nsteps
        %   compute for j+1
        g   = FEM.M*du(:,j) - FEM.Tsteps(j)*dN*HI.G(:,j+1);
        du  (FEM.p,j+1) = HI.R{j+1}\(HI.R{j+1}'\g(FEM.p,:));
    end
    
    ind     = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    Jx(:,ind)   = (FEM.C*du)*FEM.obs_int;
end

%}