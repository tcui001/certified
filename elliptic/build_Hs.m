
nlabs = 32;
parpool('local', nlabs)

par_Hg  = Composite(nlabs);
par_Hh  = Composite(nlabs);

opts.issym  = 1;
opts.isreal = 1;
myeig = @(HI) eigs(@(dv) matvec_PPGNH(model, obs, prior, HI, dv), prior.DoF,  30, 'LA', opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samples = out_dili.v_samples';
n = size(samples, 2);

% posterior

Md  = distributed(samples);
parallel_H_sum;
new_H = tmp_Hg/n;
lis_H = tmp_Hh/n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prior

M   = randn(prior.DoF, n);
Md  = distributed(M);

parallel_H_sum;
as_H = tmp_Hg/n;
prlis_H = tmp_Hh/n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Laplace

[V,d]   = eigen_PPGNH(model, obs, prior, vmap, 1E-5, 50);
s       = (d+1).^(-0.5) - 1;
r       = randn(prior.DoF, n);
M       = r + V*scale_rows(V'*r, s) + repmat(vmap, 1, n);
Md      = distributed(M);

parallel_H_sum;
lap_H = tmp_Hg/n;
laplis_H = tmp_Hh/n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_H2 = out_dili.grad2;

save H_data.mat new_H2 new_H lis_H as_H prlis_H lap_H laplis_H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
tic;
H1 = zeros(prior.DoF);
H2 = zeros(prior.DoF);

for i = 1:size(M, 2)
    [~,~,~,gmllkd,HI] = minus_log_post(model, obs, prior, M(:,i));
    
    H1 = H1 + gmllkd*gmllkd';
    
    [V, D]      = eigs(@(dv)  matvec_PPGNH(model, obs, prior, HI, dv), prior.DoF,  30, 'LA', opts);
    %[V, D]      = myeig(HI);
    d           = diag(D);
    ind         = d>1E-6;
    
    H2 = H2 + V(:,ind)*scale_cols(V(:,ind),d(ind))';
end
toc
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
M = randn(4, 100);
B = distributed(M);
Ts = Composite(nlabs);

tic
spmd
    Tl = 0;
    V = getLocalPart(B);
    
    for i = 1:size(V, 2)
        Tl = Tl+ V(:,i);
    end
    Ts = Tl;
end
toc

T1 = 0;
tic;
for i = 1:nlabs
    T1 = T1 + Ts{i};
end
toc

tic;
T2 = 0;

for i = 1:size(M, 2)
    T2 = T2 + M(:,i);
end
toc

%}