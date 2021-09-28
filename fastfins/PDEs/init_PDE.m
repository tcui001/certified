function FEM = init_PDE(model_def, obs_def, mesh)
%INIT_PDE
%
% Tiangang Cui, 11/May/2014

switch model_def.test_case
    case {'EIT'} % EIT
        FEM         = EIT_make_outputs(mesh, obs_def);
        FEM         = EIT_make_FEM(FEM, mesh);
        FEM.problem = 'Laplace';
    case {'Laplace_T1'} % test case for DILI
        FEM         = laplace_make_outputs(mesh, obs_def);
        FEM         = laplace_make_FEM_T1(FEM, mesh);
        FEM.problem = 'Laplace';
    case {'Laplace_T2'} % test case for dim redu, ex 1
        FEM         = laplace_make_outputs(mesh, obs_def);
        FEM         = laplace_make_FEM_T2(FEM, mesh);
        FEM.problem = 'Laplace';
    case {'Laplace_T3'} % test case for dim redu, ex 2
        FEM         = laplace_make_outputs(mesh, obs_def);
        FEM         = laplace_make_FEM_T3(FEM, mesh);
        FEM.problem = 'Laplace';
    case {'Laplace_TT'} % test case for dim redu, ex 2
        FEM         = laplace_make_outputs(mesh, obs_def);
        FEM         = laplace_make_FEM_TT(FEM, mesh);
        FEM.problem = 'Laplace';
    case {'Laplace_DILI'} % test case for dim redu, ex 2
        FEM         = laplace_make_outputs(mesh, obs_def);
        FEM         = laplace_make_FEM_DILI(FEM, mesh);
        FEM.problem = 'Laplace';
    case{'Heat_T1'} % test case 1
        FEM         = heat_make_outputs(mesh, model_def.time_def, obs_def);
        FEM         = heat_make_FEM_T1 (FEM,  mesh);
        FEM.problem = 'Heat';
        FEM.GMRES_flag      = model_def.use_GMRES;
    case{'Heat_T2'} % test case 1
        FEM         = heat_make_outputs(mesh, model_def.time_def, obs_def);
        FEM         = heat_make_FEM_T2 (FEM,  mesh);
        FEM.problem = 'Heat';
        FEM.GMRES_flag      = model_def.use_GMRES;
        FEM.pote_t  = model_def.param_pote_t;
        FEM.pote_b  = model_def.param_pote_b;
        FEM.fs      = FEM.ft*FEM.pote_t + FEM.fb*FEM.pote_b + FEM.f;
    case{'RD_T0'}   % test case 0
        FEM         = RD_make_outputs(mesh, obs_def);
        FEM         = RD_make_FEM(FEM,  mesh);
        FEM.problem = 'RD';
        FEM.GMRES_flag      = model_def.use_GMRES;
        FEM.fast_adj_flag   = model_def.fast_adj;
        FEM.ML_flag         = model_def.use_ML;
        FEM.a               = model_def.param_a;
        FEM.kappa           = model_def.param_k;
        FEM.dt_max          = model_def.dt_max;
        FEM.dt_init         = model_def.dt_init; 
        FEM.dt_multiplier   = model_def.dt_multiplier;
        FEM.Nmaxsteps       = 500;
end

end
