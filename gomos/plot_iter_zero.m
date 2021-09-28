
load run_iter_zero_ip.mat
pr0_dkl = zeros(n_iter, 1);
pr0_dh = zeros(n_iter, 1);
pr0_bounds = zeros(n_iter, 1);
pr0_nbasis = zeros(n_iter, 1);
for i = 1:n_iter
    [pr0_dkl(i), pr0_dh(i)] = calc_distances(full_mllkd, mllkds{i}');
    
    [V, S, ~] = svd(Hs{i});
    s = diag(S);
    pr0_nbasis(i) = prior_rs{i}.DoF;
    pr0_bounds(i) = sum(s.^2) - sum(s(1:pr0_nbasis(i)).^2);
end


load run_iter_zero.mat
pr0i_dkl = zeros(n_iter, 1);
pr0i_dh = zeros(n_iter, 1);
pr0i_bounds = zeros(n_iter, 1);
pr0i_nbasis = zeros(n_iter, 1);
for i = 1:n_iter
    [pr0i_dkl(i), pr0i_dh(i)] = calc_distances(full_mllkd, mllkds{i}');
    [V, S, ~] = svd(Hs{i});
    s = diag(S);
    pr0i_nbasis(i) = prior_rs{i}.DoF;
    pr0i_bounds(i) = sum(s.^2) - sum(s(1:pr0i_nbasis(i)).^2);
end

%{
figure('position', [100, 100, 900 350]);
subplot(1,2,1)
ax = gca;
ax.ColorOrderIndex = 3;
semilogy(1:n_iter, pr0_bounds, 'x-')
hold on
ax.ColorOrderIndex = 1;
semilogy(1:n_iter, pr0_dkl, 'v-')
ax.ColorOrderIndex = 4;
semilogy(1:n_iter, pr0i_bounds, 'o-')
ax.ColorOrderIndex = 2;
semilogy(1:n_iter, pr0i_dkl, 's-')
xlabel('iteration')
h = legend('error bound, without IP', '$D_{\rm KL}(\nu || \hat{\nu}_r^k), \; Y_i=m$, without IP', 'error bound, IP', '$D_{\rm KL}(\nu || \hat{\nu}_r^k), \; Y_i=m$, IP');
set(h, 'box', 'off', 'interpreter', 'latex');
set(gca, 'xtick', 1:5)
axis([1, 5, 1E-2 1E7])

subplot(1,2,2)
ax = gca;
ax.ColorOrderIndex = 1;
semilogy(1:n_iter, pr0_dh, 'v-')
hold on
semilogy(1:n_iter, pr0i_dh, 's-')
h = legend('$D_{\rm H}(\nu || \hat{\nu}_r^k), \; Y_i=m$, without IP', '$D_{\rm H}(\nu || \hat{\nu}_r^k), \; Y_i=m$, IP');
set(h, 'box', 'off', 'interpreter', 'latex');
xlabel('iteration')
axis([1, 5, 0.05 1])
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position', [100, 100 700 400]);
subplot(2,2,1)
semilogy(1:n_iter, pr0_bounds, 'o-')
hold on
semilogy(1:n_iter, pr0_dkl, 's-')
%yyaxis right
%plot(1:n_iter, pr0_nbasis, 'x-')
plot([1, n_iter], [1E-2, 1E-2], 'k')
h = legend('$\frac12\mathcal{R}_\Gamma(P_r, H_{ref})$', '$D_{\rm KL}(\nu || \hat{\nu}_r), \; Y_1=m$', '$10^{-2}$');
set(h, 'box', 'off', 'interpreter', 'latex');
xlabel('iteration k')
ylabel('error')
set(gca, 'xtick', 1:5)
axis([1, 5, 1E-3 1E7])
set(gca, 'ytick', 10.^(-3:2:7))
title('fix error threshold = $10^{-2}$')


subplot(2,2,3)
plot(1:n_iter, pr0_nbasis, 'o-')
xlabel('iteration k')
ylabel('rank')
set(gca, 'xtick', 1:5)
axis([1, 5, 30 40])


subplot(2,2,2)
semilogy(1:n_iter, pr0i_bounds, 'o-')
hold on
semilogy(1:n_iter, pr0i_dkl, 's-')
plot([1, n_iter], [1E-2, 1E-2], 'k')
%h = legend('error bound', '$D_{\rm KL}(\nu || \hat{\nu}_r^{(k)}), \; Y_i=m$', '$10^{-2}$');
%set(h, 'box', 'off', 'interpreter', 'latex');
set(gca, 'xtick', 1:5)
axis([1, 5, 1E-3 1E7])
set(gca, 'ytick', 10.^(-3:2:7))
title('fix rank = 30')
xlabel('iteration k')
ylabel('error')


subplot(2,2,4)
plot(1:n_iter, pr0i_nbasis, 'o-')
xlabel('iteration k')
ylabel('rank')
set(gca, 'xtick', 1:5)
axis([1, 5, 30 40])