function FEM = EIT_make_outputs(mesh, obs_def)
%EIT_MAKE_OUTPUTS   
%
% generates measurements points and indices w.r.t. the FEM mesh
%
% Tiangang Cui, 03/May/2012

tol             = 1E-10;
FEM.Nsensors    = 4*length(obs_def.locs);
FEM.sensors     = zeros(4*length(obs_def.locs),1);
count           = 0;
% y = 0
for i   = 1:length(obs_def.locs)
    j       = find(abs(mesh.node(2,:))<tol      & abs(mesh.node(1,:)-obs_def.locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

% x = xyratio
for i   = 1:length(obs_def.locs)
    j       = find(abs(mesh.node(1,:)-mesh.xyratio)<tol & abs(mesh.node(2,:)-obs_def.locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

% y = 1
for i   = length(obs_def.locs):-1:1
    j       = find(abs(mesh.node(2,:)-1)<tol    & abs(mesh.node(1,:)-obs_def.locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

% x = 0
for i   = length(obs_def.locs):-1:1
    j       = find(abs(mesh.node(1,:))<tol      & abs(mesh.node(2,:)-obs_def.locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

FEM.Ndatasets   = FEM.Nsensors;
FEM.Ndata       = FEM.Nsensors^2;
FEM.Nsteps      = 1;

FEM.C           = sparse([],[],[], FEM.Nsensors, mesh.Nnode, FEM.Nsensors); % make the observation matrix
FEM.C(:,FEM.sensors)= speye(FEM.Nsensors);

end
