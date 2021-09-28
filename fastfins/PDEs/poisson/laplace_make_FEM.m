function FEM = laplace_make_FEM(FEM, mesh)
%LAPLACE_MAKE_FEM    
%
% makes the basic  FEM structure
%
% Tiangang Cui, 03/May/2014

FEM.W1              = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.W2              = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.W3              = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.M               = sparse([],[],[],mesh.Nnode,mesh.Nnode,10*mesh.Nel);
FEM.Mb              = sparse([],[],[],mesh.Nnode,mesh.Nnode, 4*mesh.Nbndf);

for i = 1:mesh.Nbndf
    ind             = mesh.node_map_bnd(:,i);
    dx              = mesh.node(:,ind(2)) - mesh.node(:,ind(1));
    FEM.Mb(ind,ind) = FEM.Mb(ind,ind) + norm(dx)*mesh.locmass_bnd;
end

for i = 1:mesh.Nel
    ind             = mesh.node_map(:,i);
    % assume const and identity Jacobian, note here the local stiffness
    % matrix does not contain the Jacobian since square grid is used
    dx              = mesh.node(:,ind(3)) - mesh.node(:,ind(1));
    detJ            = prod(abs(dx)); % iJ = diag(1./dx);
        
    FEM.W1(ind,i)   = mesh.w1;
    FEM.W2(ind,i)   = mesh.w2;
    FEM.W3(ind,i)   = mesh.w3;
    
    FEM.M(ind,ind)  = FEM.M(ind,ind) + mesh.locmass*detJ;    
end

FEM.MR              = chol(FEM.M);

% also set up the constraint vector for esential conditions (mean on boundary is zero)
c = sparse(mesh.node_map_bnd(1,:),1,1/mesh.Nbndf,mesh.Nnode,1) + ...
    sparse(mesh.node_map_bnd(2,:),1,1/mesh.Nbndf,mesh.Nnode,1); % projection is now c*c'

FEM.c               = c*c'*mesh.Nel; % penalty matrix

% Precompute permutation vector for modified stiffness matrix to give 
% minimum bandwidth. 
Apen                = FEM.W1*FEM.W1' + FEM.W2*FEM.W2' + FEM.W3*FEM.W3' + FEM.c;

FEM.p               = symamd(Apen);   % this is a MatLab function
FEM.r(FEM.p)        = 1:size(FEM.p,2); % inverse permutation

% make the observation matrix
%FEM.C               = sparse([],[],[],FEM.Nsensors,mesh.Nnode,FEM.Nsensors);
%FEM.C(:,FEM.sensors)= speye(FEM.Nsensors);

FEM.DoF             = mesh.Nnode;
FEM.mesh            = mesh;

for i = 1:mesh.Nbndf
    ind = mesh.node_map_bnd(:,i);
    dx = mesh.node(:,ind(2)) - mesh.node(:,ind(1));
    FEM.Mb(ind,ind) = FEM.Mb(ind,ind) + norm(dx)*mesh.locmass_bnd;
end
%FEM.upscale_type = 0;

end
