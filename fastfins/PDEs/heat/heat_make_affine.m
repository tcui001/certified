function FEM = heat_make_affine(FEM, param)
%HEAT_MAKE_AFFINE  
%
% makes the basic FEM structure for affine problem
%
% Tiangang Cui, 31/Oct/2012

FEM.M_temp      = cell(param.S_DoF,1);
FEM.N_temp      = cell(param.K_DoF,1);

for i = 1:param.S_DoF
    FEM.M_temp{i}   = sparse([],[],[], FEM.mesh.Nnode, FEM.mesh.Nnode, 4*FEM.mesh.Nel);
    
    for j = 1:FEM.mesh.Nel
        ind         = FEM.mesh.node_map(:,j);
        dx          = mesh.node(:,ind(3)) - mesh.node(:,ind(1));
        detJ        = prod(abs(dx));
        FEM.M_temp{i}(ind,ind)  = FEM.M_temp{i}(ind,ind) + FEM.mesh.locmass*detJ*param.S_basis(j,i);
    end
end

for i = 1:param.K_DoF
    FEM.N_temp{i}   = sparse([],[],[], FEM.mesh.Nnode, FEM.mesh.Nnode, 4*FEM.mesh.Nel);
    
    for j = 1:FEM.mesh.Nel
        ind         = FEM.mesh.node_map(:,j);
        FEM.N_temp{i}(ind,ind)  = FEM.N_temp{i}(ind,ind) + FEM.mesh.locstiff*param.K_basis(j,i);
    end
end

end