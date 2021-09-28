function FEM = RD_make_outputs(mesh, obs_def)
%HEAT_MAKE_OUTPUTS
%
% generates measurements points and indices w.r.t. the FEM mesh
%
% Tiangang Cui, 03/May/2014

tol             = 1E-10;

switch  obs_def.type
    case {1}
        [xx,yy] = meshgrid(obs_def.locs, obs_def.locs);
        locs    = [xx(:) yy(:)];
        %ind = rand(size(temp,1),1)<1.1;
        %locs = temp(ind,:);
        
        FEM.sensors         = zeros(size(locs,1),1);
        
        for i   = 1:size(locs,1)
            ind = find(sum(abs(mesh.node - locs(i,:)'*ones(1,mesh.N_node)))<tol);
            FEM.sensors(i)  = ind;
        end
        
    case {2}
        ind     = mesh.node(1,:)>(obs_def.locs(1)+tol) & mesh.node(1,:)<(obs_def.locs(2)-tol) & ...
                  mesh.node(2,:)>(obs_def.locs(3)+tol) & mesh.node(2,:)<(obs_def.locs(4)-tol);
        FEM.sensors         = find(ind);
        
    case {3}
        ind     = mesh.node(1,:)>tol & mesh.node(1,:)<(1-tol) & ...
                  mesh.node(2,:)>tol & mesh.node(2,:)<(1-tol);
        FEM.sensors         = find(ind);
end

FEM.Nsensors    = length(FEM.sensors);
FEM.Tfinal      = obs_def.Tfinal;
FEM.Ndatasets   = obs_def.Ntime;


end