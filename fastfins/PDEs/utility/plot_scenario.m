% plot the true conductivity, observation sites, the true solution, and the
% data
% Tiangang Cui, 21/May/2012

switch model_def.problem
    case {'Log_Gaussian_Cox'}
        figure('position',[100 100 700 300])
        subplot(1,2,1)
        meshc(mesh_dual.gx,mesh_dual.gy,reshape(param.true_x*FEM.a,mesh.N_side,mesh.N_side))
        title('True Conductivity')
        xlabel('x')
        ylabel('y')
        
        % plot the true solution
        subplot(1,2,2)
        meshc(mesh_dual.gx,mesh_dual.gx,reshape(obs.d,mesh.N_side,mesh.N_side));
        title('Data')
        xlabel('x')
        ylabel('y')
    case {'Laplace'}
        figure('position',[100 100 700 600])
        subplot(2,2,1)
        meshc(mesh_dual.gx,mesh_dual.gy,reshape(param.true_x,mesh.N_side,mesh.N_side))
        title('True Conductivity')
        xlabel('x')
        ylabel('y')
        
        % plot the true solution
        subplot(2,2,2)
        meshc(mesh.gx,mesh.gx,reshape(soln.G(:,1),mesh.N_side+1,mesh.N_side+1));
        title('True Solution')
        xlabel('x')
        ylabel('y')
        
        % plot the observation sites
        subplot(2,2,3)
        plot_grid(mesh,'b',':');
        hold on
        plot(mesh.node(1,obs.elect),mesh.node(2,obs.elect),'bo')
        title('Observation Site')
        xlabel('x')
        ylabel('y')
        
        % plot the observation data
        subplot(2,2,4)
        plot(1:obs.N_obs,obs.d(:),'ro',1:obs.N_obs,soln.d(:),'bx-');
        title('Data')
        xlabel('x')
        ylabel('y')
    case {'Heat_single'}
        figure('position',[100 100 700 600])
        subplot(2,2,1)
        meshc(mesh_dual.gx,mesh_dual.gy,reshape(param.true_x,mesh.N_side,mesh.N_side))
        title('True Conductivity')
        xlabel('x')
        ylabel('y')
        
        % plot the true solution
        subplot(2,2,2)
        meshc(mesh.gx,mesh.gx,reshape(soln.G(:,end),mesh.N_side+1,mesh.N_side+1));
        title('True Solution')
        xlabel('x')
        ylabel('y')
        
        % plot the observation sites
        subplot(2,2,3)
        plot_grid(mesh,'b',':');
        hold on
        plot(mesh.node(1,obs.elect),mesh.node(2,obs.elect),'bo')
        title('Observation Site')
        xlabel('x')
        ylabel('y')
        
        % plot the observation data
        subplot(2,2,4)
        plot(obs.ts,obs.d','ro',obs.ts,soln.d','bx-');
        title('Data')
        xlabel('x')
        ylabel('y')
    case {'Heat_double'}
        figure('position',[100 100 700 600])
        subplot(2,2,1)
        meshc(mesh_dual.gx,mesh_dual.gy,reshape(param.true_x(param.K_ind),mesh.N_side,mesh.N_side))
        title('True Peameability')
        xlabel('x')
        ylabel('y')
        
        % plot the true solution
        subplot(2,2,2)
        meshc(mesh_dual.gx,mesh_dual.gy,reshape(param.true_x(param.S_ind),mesh.N_side,mesh.N_side))
        title('True Storage')
        xlabel('x')
        ylabel('y')
        
        % plot the observation sites
        subplot(2,2,3)
        plot_grid(mesh,'b',':');
        hold on
        plot(mesh.node(1,obs.elect),mesh.node(2,obs.elect),'bo')
        title('Observation Site')
        xlabel('x')
        ylabel('y')
        
        % plot the observation data
        subplot(2,2,4)
        plot(obs.ts,obs.d,'ro',obs.ts,soln.d,'bx-');
        title('Data')
        xlabel('x')
        ylabel('y')
end