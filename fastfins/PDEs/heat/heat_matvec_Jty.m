function Jty = heat_matvec_Jty(FEM, HI, dy)
%HEAT_MATVEC_JTY
%
% Adjoint solve
%
% Tiangang Cui, 09/May/2014 

TOL         = 1E-6;
MAXIT       = 50;

N           = size(dy, 2)/FEM.Ndatasets;
Jty         = zeros( HI.NoP, N);
lambda      = zeros(FEM.DoF, 1);

for i   = 1:N
    ind         = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    dU          = FEM.C'*(dy(:,ind)*FEM.obs_int');   
    lambda(:)   = 0;
    
    for j = (FEM.Nsteps+1):-1:2
        g       = FEM.M*lambda - dU(:,j);
        
        if FEM.GMRES_flag
            B   = FEM.M + FEM.Tsteps(j-1)*HI.N;
            [lambda,~]      = minres(B, g, TOL, MAXIT, HI.L{j}, HI.L{j}');
        else
            lambda(FEM.p)   = HI.L{j}'\(HI.L{j}\g(FEM.p));
        end
        
        T1      = reshape(lambda(FEM.mesh.node_map),4,FEM.mesh.Nel);
        T2      = reshape(HI.G(FEM.mesh.node_map,j),4,FEM.mesh.Nel);
        
        Jty(:,i)= Jty(:,i) + sum(T1.*(FEM.mesh.locstiff*T2))'*FEM.Tsteps(j-1);
    end

end

end

%{

    for j = (FEM.Nsteps+1):-1:2
        if  j   == (FEM.Nsteps+1)
            g   = -dU(:,j);
        else
            g   = FEM.M*lambda(:,j) - dU(:,j);
        end
        
        if FEM.GMRES_flag
            B   = FEM.M + FEM.Tsteps(j-1)*FEM.N;
            lambda(:,j)     = minres(B, g, TOL, MAXIT, soln.L{i+1}, soln.L{i+1}');
        else
            lambda(FEM.p,j) = HI.L{j}'\(HI.L{j}\g);
        end
    end
    
    % step one
    lambda(:,1) = FEM.M*lambda(:,2) - dU(:,1);
    
    for j = 2:(FEM.Nsteps+1)
        T1      = reshape(lambda(FEM.mesh.node_map,j),4,FEM.mesh.Nel);
        T2      = reshape(  HI.G(FEM.mesh.node_map,j),4,FEM.mesh.Nel);
        Jty(:,i)= Jty(:,i) + sum(T1.*(FEM.mesh.locstiff*T2))'*FEM.Tsteps(j-1);
    end

%}


%{

TOL         = 1E-6;
MAXIT       = 50;

N           = size(dy, 2)/FEM.Ndatasets;
Jty         = zeros(FEM.DoF, N);
%lambda     = zeros(FEM.DoF*2, HI.Nend+1); no storage for lambda
dU          = zeros(FEM.DoF*2, HI.Nend+1);
MLJlambda   = zeros(FEM.DoF*2, 1);

% form the dResdalpha as a vector
tmp_u       = HI.G(FEM.ind_u,:);
dResdalpha  = scale_cols(tmp_u.*(1-tmp_u), HI.dts(1:HI.Nend+1));


for i   = 1:N
    % perturbation on the data space
    ind             = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    dU(:,:)         = 0;
    dU(FEM.ind_u,:) = FEM.C'*(dy(:,ind)*HI.T');    
    
    MLJlambda(:)    = 0;
    
    % adjoint    
    for j = (HI.Nend+1):-1:2
        dt          = HI.dts(j);
        mudt        = FEM.mu*dt;
        epdt        = FEM.epsilon*dt;
        adt         = FEM.a*dt;
        bdt         = FEM.b*dt;

        A           = kron(diag([mudt, epdt]), FEM.iMLK) + kron([1, dt; -adt, 1+bdt], FEM.I);
        J           = direct_dResds(FEM, HI.alpha, A, dt, HI.G(:,j));
        
        g           = MLJlambda + dU(:,j); % RHS
        if HI.pre_flag
            if FEM.GMRES_flag
                H       = J'*FEM.ML2*J; 
                lambda  = minres(H, g, TOL, MAXIT, HI.Ls{j}, HI.Ls{j}');
            else
                lambda  = HI.Ls{j}'\(HI.Ls{j}\g);
            end
        else
            H           = J'*FEM.ML2*J;
            if FEM.GMRES_flag
                L       = make_precond_sym(H);
                lambda  = minres(H, g, TOL, MAXIT, L, L');
            else
                lambda  = H\g;
            end
        end
        
        MLJlambda   = FEM.vecML.*(J*lambda);
        Jty(:,i)    = Jty(:,i) + dResdalpha(:,j).*MLJlambda(FEM.ind_u);
    end
end

%}