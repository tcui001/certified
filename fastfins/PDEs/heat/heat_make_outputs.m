function FEM = heat_make_outputs(mesh, time_def, obs_def)
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
    case {4}
        FEM.sensors = zeros(size(obs_def.locs,1),1);
        for i = 1:size(obs_def.locs,1)
            ind     = find(sum(abs(mesh.node - obs_def.locs(i,:)'*ones(1,mesh.Nnode)))<tol);
            FEM.sensors(i)  = ind;
        end
    case {5}
        N           = size(obs_def.locs,1);
        FEM.C       = [];
        for k = 1:N
            ind     = mesh.centers(1,:)>(obs_def.locs(k,1)+tol) & mesh.centers(1,:)<(obs_def.locs(k,2)-tol) & ...
                      mesh.centers(2,:)>(obs_def.locs(k,3)+tol) & mesh.centers(2,:)<(obs_def.locs(k,4)-tol);
            tmp     = 1:mesh.Nel;
            elems   = tmp(ind);
            C       = sparse([],[],[],1,mesh.Nnode,ceil((sqrt(sum(ind))+1)^2));
            for i = 1:sum(ind)
                nodes   = mesh.node_map(:,elems(i));
                dx      = mesh.node(:,nodes(3)) - mesh.node(:,nodes(1));
                detJ    = prod(abs(dx)); % iJ = diag(1./dx);
                locs    = 0.25*detJ*ones(1,4);
                C(1,nodes)  = C(1,nodes) + locs;
            end
            FEM.C   = [FEM.C; C];
        end
        
end

if  obs_def.type == 5
    FEM.Nsensors        = size(FEM.C, 1);
else
    FEM.Nsensors        = length(FEM.sensors);
    FEM.C               = sparse([],[],[], FEM.Nsensors, mesh.Nnode, FEM.Nsensors); % make the observation matrix
    FEM.C(:,FEM.sensors)= speye(FEM.Nsensors);
end
FEM.Tstart              = obs_def.Tstart;
FEM.Tfinal              = obs_def.Tfinal;
FEM.Ndatasets           = obs_def.Ntime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FEM.Nobstime    = obs_def.Ntime;
%FEM.Tsteps      = linspace(0, obs_def.Tfinal, obs_def.Ntime);

% process the time steps for the FEM, use time_def
%Tsteps          = [time_def.Tinit*ones(1, time_def.NTinit), ...
%                   logspace(log10(time_def.Tinit), log10(time_def.Tnormal), time_def.NTtransit), ...
%                   time_def.Tnormal*ones(1, time_def.NTnormal)]; % setup the time
               
Tsteps          = [logspace(log10(time_def.Tinit), log10(time_def.Tnormal), time_def.NTtransit), ...
                   time_def.Tnormal*ones(1, time_def.NTnormal)]; % setup the time
                   
FEM.Tsteps      = Tsteps/sum(Tsteps)*FEM.Tfinal;
FEM.Tsteps(end) = FEM.Tsteps(end)+time_def.Tinit;
FEM.Nsteps      = length(FEM.Tsteps);

% process observations 
obs_Tsteps      = linspace(0, FEM.Tfinal, FEM.Ndatasets);
FEM.obs_int     = sparse([],[],[], FEM.Nsteps, FEM.Ndatasets, FEM.Ndatasets*2);
t_start         = 0;
for           i = 1:FEM.Nsteps
    t_end       = t_start + FEM.Tsteps(i);
    t_ind       = find(obs_Tsteps>t_start & obs_Tsteps<=t_end);
    temp1       = (t_end - obs_Tsteps(t_ind))/FEM.Tsteps(i);   % weighting at t_start
    temp2       = (obs_Tsteps(t_ind) - t_start)/FEM.Tsteps(i); % weighting at t_end
    t_start     = t_end; % increment the time
    
    FEM.obs_int(i,  t_ind)  = temp1;
    FEM.obs_int(i+1,t_ind)  = temp2;    
end

end