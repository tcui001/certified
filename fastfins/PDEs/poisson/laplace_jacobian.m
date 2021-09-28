function J = laplace_jacobian(FEM, soln)
%
% compute the Jacobian
% Tiangang Cui, 10/Mar/2013

dual(FEM.p,:)   = -full(soln.upper\(soln.lower\FEM.C(:,FEM.p)'));
J               = zeros(FEM.Nsensors*FEM.Ndatasets, FEM.mesh.Nel);

for j = 1:FEM.Ndatasets
    bi          = (j-1)*FEM.Nsensors;
    T2          = FEM.mesh.locstiff*reshape(soln.G(FEM.mesh.node_map,j), 4, FEM.mesh.Nel);
    for i = 1:FEM.Nsensors
        T1          = reshape(dual(FEM.mesh.node_map,i), 4, FEM.mesh.Nel);
        J(bi+i,:)   = J(bi+i,:) + sum(T1.*T2);
    end
end
