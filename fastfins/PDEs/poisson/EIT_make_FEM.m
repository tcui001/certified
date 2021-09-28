function FEM = EIT_make_FEM(FEM, mesh)
%EIT_MAKE_FEM      
%
% makes the FEM structure for the EIT problem
%
% Tiangang Cui, 03/May/2014

FEM     = laplace_make_FEM(FEM, mesh);
FEM.fs  = sparse([],[],[], mesh.Nnode, FEM.Nsensors, mesh.Nnode*FEM.Nsensors); 

for i   = 1:FEM.Nsensors
    % implement boundary reference, i.e. current removed evenly round
    % boundary
    % FEM.fs(mesh.node_map_bnd(1,:),count) = -(1/mesh.N_bnd_f)*ones(mesh.N_bnd_f,1); 
    % fsp(FEM.sensors(count),count)        = FEM.fs(FEM.sensors(count),count) + 1;
    
    fsm = zeros(mesh.Nnode, 1);
    
    fsm(mesh.node_map_bnd(1,:)) = -(1/4)*ones(mesh.Nbndf,1); 
    fsm = FEM.Mb*fsm;
    
    xn  = mesh.node(:,FEM.sensors(i));
    dis = sum( (mesh.node - repmat(xn, 1, mesh.Nnode)).^2 )';
    fsp = FEM.M*exp(-0.5*dis/(2*0.02^2))/(pi*0.02^2);
    
    %sum(fsp)
    %sum(fsm)
    
    FEM.fs(:,i)  = fsp + fsm;
end

end