function [model, obs, prior] = setup_PDE(model_def, obs_def, prior_def, output)
%SETUP_PDE
%
% Tiangang Cui, 11/May/2014

gx              = linspace(0, model_def.xyratio,  model_def.xyratio*model_def.mesh_size+1);
gy              = linspace(0, 1,                  model_def.mesh_size+1);
[fmesh, pmesh]  = rectmesh2d(gx,gy);  % mesh maker
fmesh.xyratio   = model_def.xyratio;
pmesh.xyratio   = model_def.xyratio;

model           = init_PDE(model_def, obs_def, fmesh); % setup model, and dimension of the observations
prior           = init_prior_dist(prior_def, pmesh);   % setup prior

%%%%%%%%%%%%%%%%%%% loading test parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prior.true_x    = load_test_image(prior, prior_def.true_image, output);

%%%%%%%%%%%%%%%%%%% end of loading test param %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.soln      = forward_solve(model, prior.true_x);  % reference solution

%%%%%%%%%%%%%%%%%%% loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  exist(output.data_str,'file')
    obs_ref     = load(output.data_str);
    obs.data    = obs_ref.data;
    obs.std     = obs_ref.std;
else
    % the s.t.d. is calculated from the signal to noise ratio
    if obs_def.s2n > 0
        obs.std = max(abs(model.soln.d(:)))/obs_def.s2n;
    else
        obs.std = obs_def.std;
    end
    % generate data
    obs.data    = model.soln.d + randn(model.Nsensors, model.Ndatasets)*obs.std;
    save(output.data_str, '-struct', 'obs');
end

obs.Nsensors    = model.Nsensors;
obs.Ndatasets   = model.Ndatasets;
obs.Ndata       = obs.Nsensors*obs.Ndatasets;

%%%%%%%%%%%%%%%%%%% end of loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end