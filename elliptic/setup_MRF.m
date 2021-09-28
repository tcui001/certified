% setup the model problem, with Markov random field prior
% Tiangang Cui 17/May/2014

%**************************************************************************
% options for setup the PDE model  
%**************************************************************************

model_def.test_case     = 'Laplace_TT'; % test case
model_def.mesh_size     = 40;           % mesh size on y dimension
model_def.xyratio       = 3;            % xyratio of the mesh
model_def.use_GMRES     = false;

prior_def.func          = 'log'; % also can use 'erf', if not given, default in no transformation
prior_def.log_thres     = 0;     % this must be set for log transformation
prior_def.cov_type      = 'MRF'; % MRF setup
[Q,~]   = qr([1 1; 1 -1]);
corr    = Q*diag([10, 1]/2)*Q';
prior_def.cond          = [corr(1,1), corr(2,2), corr(1,2)];
prior_def.k             = 2;
prior_def.sigma         = 2;
prior_def.mean          = log(1E-3);

prior_def.true_image.type   = 'Prior'; % also can choose 'CF' for EIT and DILI
prior_def.true_image.base   = 2;
prior_def.true_image.range  = 2;

%**************************************************************************
% options for setup the prior and the ``true'' parameter 
%**************************************************************************

obs_def.type            = 4;
dx                      = 0.2;
[x1, y1]                = meshgrid(0.3:dx:0.7, 0.2:dx:0.7);
[x2, y2]                = meshgrid(2.4:dx:2.6, 0.4:dx:0.6);
obs_def.locs            = [x1(:) y1(:); x2(:) y2(:)];
obs_def.s2n             = 20; 

%**************************************************************************
% output options
%**************************************************************************

output.data_str         = [pwd '/MRF_snr20_data_new.mat'];
output.image_str        = [pwd '/MRF_image_new.mat'];

%**************************************************************************
% setup model
%**************************************************************************
[model, obs, prior] = setup_PDE(model_def, obs_def, prior_def, output);

vmap                    = get_map_matlab(model, obs, prior, randn(prior.DoF, 1)); 

umap    = matvec_prior_L(prior, vmap) + prior.mean_u;
xmap    = exp(umap);
tic; 
model.GMRES_flag = 0;
soln    = forward_solve(model, xmap);
toc

tic;
[V, d] = eigen_PPGNH(model, obs, prior, vmap, 1E-6, 50);
toc