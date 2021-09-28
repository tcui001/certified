% debug the grad and Hessmult
forward = FEM;
prior.beta = 0;
% generate a test parameter from the prior
z = dist_random(prior,1);

% param transformation
%[x dxdy dx2dy2] = param_z2x( param, z );

% forward model run
[f g hessinfo] = minus_log_post(mesh, obs, prior, param, forward, z);

tol = 1e-6;
% finite difference gradient
f1 = zeros(size(z));
f2 = zeros(size(z));
for i = 1:param.DoF_GP
    
    perturb = zeros(size(z));
    perturb(i) = tol;
    
    z1 = z-perturb;
    f1(i) = minus_log_post_simple(obs, prior, param, forward, z1);
    
    z2 = z+perturb;
    f2(i) = minus_log_post_simple(obs, prior, param, forward, z2);
    
end

g_fd = (f2-f1)/(2*tol);
figure('position',[100 100 700 300])
subplot(1,3,1)
plot(g,g_fd)
axis tight
subplot(1,3,2)
semilogy(abs(g-g_fd))
subplot(1,3,3)
semilogy(abs((g-g_fd)./g))


% finite difference Jacobian
d1 = zeros(obs.Nt*obs.N_obs,param.DoF_GP);
d2 = zeros(obs.Nt*obs.N_obs,param.DoF_GP);
fd = zeros(mesh.N_el,2);

for i = 1:param.DoF_GP
    
    perturb = zeros(size(z));
    perturb(i) = tol;
    
    z1 = z-perturb;
    x1 = param_z2x_simple(param, z1);
    soln1 = forward_solve(forward, x1);
    d1(:,i) = soln1.d(:);
    
    z2 = z+perturb;
    x2 = param_z2x_simple(param, z2);
    soln2 = forward_solve(forward, x2);
    d2(:,i) = soln2.d(:);
    
end

x = param_z2x_simple(param, z);
soln = forward_solve(forward, x);

J_fd = (d2-d1)/(2*tol);
GN_fd = J_fd'*J_fd;
g_GN_fd = J_fd'*(soln.d(:)-obs.d(:));

% Hessmult GN
forward.hess_flag = 'GN';
for i = 1:10
    dz = rand(size(z))*0.4-0.2;
    w = minus_log_post_hessmult(mesh, obs, prior, param, forward, hessinfo, dz);
    
    % finite difference Hessmult
    ndz = norm(dz);
    ddz = dz/ndz;
    
    z1 = z-tol*ddz;
    [f_temp g1] = minus_log_post(mesh, obs, prior, param, forward, z1);
    z2 = z+tol*ddz;
    [f_temp g2] = minus_log_post(mesh, obs, prior, param, forward, z2);
    
    w_fd = (g2-g1)/(2*tol)*ndz;
    
    % test the Hessmult
    % dx = dxdy*dz;
    wz = GN_fd * dz;
    w_GN_fd  = wz + prior.beta*prior.P*dz;
    
    figure('position',[100 100 700 300])
    subplot(1,3,1)
    plot(w,w_GN_fd)
    axis tight
    subplot(1,3,2)
    semilogy(abs(w-w_GN_fd))
    subplot(1,3,3)
    semilogy(abs((w-w_GN_fd)./w))
end