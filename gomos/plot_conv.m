% plot_conv

% new_* is the new method, gradient over posterior
% lis_* is the LIS over posterior
% as_* is the gradient over prior
% prlis_* is the LIS over prior
% lap_* is the gradient over laplace
% laplis_* is the LIS over laplace

% load data for kl and hellinger distance
%load mllkd_data.mat

summarise_distances

%%%%%%%%%%%%

figure('position', [100 100 900 600]);
subplot(2,4,1)
semilogy(nps, new_bounds, 'o-');
hold on
semilogy(nps, new_dklf, 'x-');
semilogy(nps, new_dkle, 's-');
semilogy(nps, new_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('$\nabla \log f, \;\rho=\nu$')
xlabel('rank of projector')
h = legend('$\frac12\mathcal{R}_\Gamma(P_r^{ref}, H_{ref})$', ...
    '$D_{\rm KL}(\nu || \hat{\nu}_r), \; Y_1=m$',...
    '$D_{\rm KL}(\nu || \hat{\nu}_r), \; M=10$',...
    '$D_{\rm KL}(\nu || \hat{\nu}_r), \; M=1$');
set(h, 'box', 'off', 'interpreter', 'latex', 'location', 'southwest');
h = get(gca,'Children');
set(gca,'Children',flipud(h))

subplot(2,4,2)
semilogy(nps, lap_bounds, 'o-');
hold on
semilogy(nps, lap_dklf, 'x-');
semilogy(nps, lap_dkle, 's-');
semilogy(nps, lap_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('$\nabla \log f, \;\rho={\rm Laplace}(\nu)$')
h = get(gca,'Children');
set(gca,'Children',flipud(h))
xlabel('rank of projector')

subplot(2,4,3)
semilogy(nps, as_bounds, 'o-');
hold on
semilogy(nps, as_dklf, 'x-');
semilogy(nps, as_dkle, 's-');
semilogy(nps, as_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('$\nabla \log f, \;\rho=\mu$')
h = get(gca,'Children');
set(gca,'Children',flipud(h))
xlabel('rank of projector')

subplot(2,4,5)
semilogy(nps, lis_bounds, 'o-');
hold on
semilogy(nps, lis_dklf, 'x-');
semilogy(nps, lis_dkle, 's-');
semilogy(nps, lis_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('$\nabla G, \;\rho=\nu$')
ylabel('KL divergence')
xlabel('rank of the projector')
h = get(gca,'Children');
set(gca,'Children',flipud(h))

subplot(2,4,6)
semilogy(nps, laplis_bounds, 'o-');
hold on
semilogy(nps, laplis_dklf, 'x-');
semilogy(nps, laplis_dkle, 's-');
semilogy(nps, laplis_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('$\nabla G, \;\rho={\rm Laplace}(\nu)$')
h = get(gca,'Children');
set(gca,'Children',flipud(h))
xlabel('rank of projector')

subplot(2,4,7)
semilogy(nps, prlis_bounds, 'o-');
hold on
semilogy(nps, prlis_dklf, 'x-');
semilogy(nps, prlis_dkle, 's-');
semilogy(nps, prlis_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('$\nabla G, \;\rho=\mu$')
h = get(gca,'Children');
set(gca,'Children',flipud(h))
xlabel('rank of projector')

subplot(2,4,4)
semilogy(nps, pr_bounds, 'o-');
hold on
semilogy(nps, pr_dklf, 'x-');
semilogy(nps, pr_dkle, 's-');
semilogy(nps, pr_dkl, 'Color', [0.8 0.8 0.8]);
axis([10, 50, 1E-2 1E4])
title('prior-based')
h = get(gca,'Children');
set(gca,'Children',flipud(h))


%%%%%%%%%%%%

figure('position',  [100 100 900 350])
subplot(1,2,1);
semilogy(nps, new_dkle, '.-', nps, lap_dkle, 'x-', nps, as_dkle, '^-');
h1 = legend('$\nabla \log f, \;\rho=\nu$', '$\nabla \log f, \;\rho={\rm Laplace}(\nu)$', '$\nabla \log f, \;\rho=\mu$'); %display legend 1
hold on
semilogy(nps, lis_dkle, 'o-', nps, laplis_dkle, 's-', nps, prlis_dkle, 'v-', nps, pr_dkle, '--');
axis([10 50 0.05 1E3])
%ylabel('KL divergence')
h1.String = {'$\nabla \log f, \;\rho=\nu$', '$\nabla \log f, \;\rho={\rm Laplace}(\nu)$', '$\nabla \log f, \;\rho=\mu$'};
set(h1, 'box', 'off', 'interpreter', 'latex', 'location', 'southwest');
xlabel('rank of projectors')
title('$D_{\rm KL}(\nu || \hat{\nu}_r), \; M=10$')

subplot(1,2,2)
p1 = semilogy(nps, new_dklf, '.-', nps, lap_dklf, 'x-', nps, as_dklf, '^-');
hold on
p2 = semilogy(nps, lis_dklf, 'o-', nps, laplis_dklf, 's-', nps, prlis_dklf, 'v-', nps, pr_dklf, '--');
h2 = legend(p2, '$\nabla G, \;\rho=\nu$', '$\nabla G, \;\rho={\rm Laplace}(\nu)$', '$\nabla G, \;\rho=\mu$', 'prior-based');            %display legend 2
set(h2, 'box', 'off', 'interpreter', 'latex', 'location', 'southwest');
axis([10 50 0.05 1E3])
title('$D_{\rm KL}(\nu || \hat{\nu}_r), \; Y_1=m$')
xlabel('rank of projectors')

%%%%%%%%%%%%

figure('position',  [100 100 900 350])
subplot(1,2,1)
semilogy(nps, new_dhe, '.-');
hold on
semilogy(nps, lap_dhe, 'x-');
semilogy(nps, as_dhe, '^-');
semilogy(nps, lis_dhe, 'o-');
semilogy(nps, laplis_dhe, 's-');
semilogy(nps, prlis_dhe, 'v-');
semilogy(nps, pr_dhe, '--');
title('$D_{\rm H}(\nu || \hat{\nu}_r), \; M=10$')
xlabel('rank of projectors')
h = legend('$\nabla \log f, \;\rho=\nu$', '$\nabla \log f, \;\rho={\rm Laplace}(\nu)$', '$\nabla \log f, \;\rho=\mu$', ...
    '$\nabla G, \;\rho=\nu$', '$\nabla G, \;\rho={\rm Laplace}(\nu)$', '$\nabla G, \;\rho=\mu$', 'prior-based');
set(h, 'box', 'off', 'interpreter', 'latex', 'location', 'southwest');
axis([10 50 1E-1 1])
set(gca, 'ytick', 0:0.1:1)

subplot(1,2,2)
semilogy(nps, new_dhf(1,:), '.-');
hold on
semilogy(nps, lap_dhf(1,:), 'x-');
semilogy(nps, as_dhf(1,:), '^-');
semilogy(nps, lis_dhf(1,:), 'o-');
semilogy(nps, laplis_dhf(1,:), 's-');
semilogy(nps, prlis_dhf(1,:), 'v-');
semilogy(nps, pr_dhf(1,:), '--');
axis([10 50 1E-1 1])
set(gca, 'ytick', 0:0.1:1)
title('$D_{\rm H}(\nu || \hat{\nu}_r), \; Y_1=m$')
xlabel('rank of projectors')