function Jx = RD_matvec_Jx(FEM, HI, dx)
%HEAT_MATVEC_JX
%
% Adjoint solve
%
% Tiangang Cui, 09/May/2014 

N   = size(dx, 2);
du  = zeros(FEM.DoF,      FEM.NSteps+1);
Jx  = zeros(FEM.Nsensors, FEM.Ndatasets*N);

for i = 1:N
    dK      = sparse(1:size(dx,1),1:size(x,1),dx(:,i));
    dN      = FEM.W1*dK*FEM.W1' + FEM.W2*dK*FEM.W2' + FEM.W3*dK*FEM.W3';
    
    du(:,1) = zeros(FEM.DoF,1);
    for j   = 1:FEM.NSteps
        %   compute for j+1
        g   = hessinfo.M*du(:,j) - FEM.Tsteps(j)*dN*HI.G(:,j+1);
        du  (FEM.p,j+1) = HI.R{j+1}\(HI.R{j+1}'\g(FEM.p,:));
    end
    
    ind     = (i-1)*FEM.Ndatasets + (1:FEM.Ndatasets);
    Jx(:,ind)   = (FEM.C*du)*FEM.obs_int;
end

end