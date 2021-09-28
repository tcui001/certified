function Jty = laplace_matvec_Jty(FEM, hessinfo, dy)
% LAPLACE_FOM_ADJOINT_JACMULT_LEFT left multiplication of the Jacobian by
%                                  perturb the data
%
%%%%%%%%%%%%%%%%%%%% Structure of the Hessinfo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tiangang Cui, 03/May/2012


%global adj_count

D               = FEM.C'*dy;  % for multiple RHS, size(dy, 2) = N * N_RHS
lambda(FEM.p,:) = -hessinfo.upper\(hessinfo.lower\D(FEM.p,:));

%adj_count = adj_count + size(D,2);

N               = size(dy,2)/FEM.Ndatasets;
Jty             = zeros(FEM.mesh.Nel, N);

for k = 1:N
    bi          = (k-1)*FEM.Ndatasets;
    for j       = 1:FEM.Ndatasets
        T1      = reshape(lambda(FEM.mesh.node_map,bi+j) ,4,FEM.mesh.Nel);
        T2      = reshape(hessinfo.G(FEM.mesh.node_map,j),4,FEM.mesh.Nel);
        Jty(:,k)    = Jty(:,k) + sum(T1.*(FEM.mesh.locstiff*T2))';
    end
end

end


