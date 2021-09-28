function vmap = get_map_matlab(model, obs, prior, v_init)
%GET_MAP_MATLAB   
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, 04/July/2020

opt  = optimoptions('fminunc', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient',true,... 
    'HessianMultiplyFcn',@(HI,v) matvec_hessian(model, obs, prior, HI, v), ...
    'Display','iter', 'MaxIterations', 500);

vmap   = fminunc_2020a(@(v) obj(model, obs, prior, v), v_init, opt);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g, hessinfo] = obj(model, obs, prior, v)
%MINUS_LOG_POST_MAP
% 
% Compute ths log posterior for optimization
%
% Tiangang Cui, 23/Mar/2013

[f, dummy, g, dummy, hessinfo] = minus_log_post(model, obs, prior, v);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = matvec_hessian(model, obs, prior, HI, dv)

w = matvec_PPGNH(model, obs, prior, HI, dv) + dv;

end
