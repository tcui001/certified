
[V,d]   = eigen_PPGNH(model, obs, prior, vmap, 1E-5, 50);
s       = (d+1).^(-0.5) - 1;
r       = randn(prior.DoF, out_dili.size);
M       = r + V*scale_rows(V'*r, s) + repmat(vmap, 1, out_dili.size);
Md      = distributed(M);

parallel_H_sum;
lap_H = tmp_Hg/out_dili.size;
laplis_H = tmp_Hh/out_dili.size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_H2 = out_dili.grad2;

save H_data2.mat new_H2 new_H lis_H as_H prlis_H lap_H laplis_H

build_redu_params;

run_reduced;

