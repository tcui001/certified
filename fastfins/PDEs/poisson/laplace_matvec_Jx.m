function Jx = laplace_matvec_Jx(FEM, HI, dx)
% LAPLACE_FOM_ADJOINT_JACMULT_RIGHT  right multiplication of the Jacobian by
%                                    perturb the data
%
%%%%%%%%%%%%%%%%%%%% Structure of the HI %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tiangang Cui, 03/May/2012

% global forward_count

N   = size(dx,2);
Jx  = zeros(FEM.Nsensors, FEM.Ndatasets*N);
du  = zeros(FEM.DoF,      FEM.Ndatasets);

for i = 1:N
    % forward_count = forward_count + 1;
    DX              = sparse(1:size(dx,1),1:size(dx,1),dx(:,i));
    dA              = FEM.W1*DX*FEM.W1' + FEM.W2*DX*FEM.W2' + FEM.W3*DX*FEM.W3';
    T1              = dA*HI.G;    % solve du from A(dx)u + Adu = 0
    du(FEM.p,:)     = -HI.upper\(HI.lower\T1(FEM.p,:));
    ind             = (1:FEM.Ndatasets) + (i-1)*FEM.Ndatasets;
    Jx(:,ind)       = FEM.C*du;
end
