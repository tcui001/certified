
%load run_iter_zero_ip.mat
pr0_dkl = zeros(n_trial, n_iter);
pr0_dh = zeros(n_trial, n_iter);
pr0_nbasis = zeros(n_trial, n_iter);
for kk = 1:n_trial
    for i = 1:n_iter
        [pr0_dkl(kk, i), pr0_dh(kk, i)] = calc_distances(full_mllkd, mllkds{kk,i}');
        pr0_nbasis(kk, i) = prior_rs{kk, i}.DoF;
    end
end


figure;

subplot(2,1,1)
plot(1:n_iter, pr0_nbasis, 'o-')
xlabel('iteration')
ylabel('basis dimension')
subplot(2,1,2)
semilogy(1:n_iter, pr0_dkl, 'o-')
xlabel('iteration')
ylabel('KL divergence')

%{

load run_iter_lap.mat
lap0_dkl = zeros(n_trial, n_iter);
lap0_dh = zeros(n_trial, n_iter);
lap0_nbasis = zeros(n_trial, n_iter);
for kk = 1:n_trial
    for i = 1:n_iter
        [lap0_dkl(kk, i), lap0_dh(kk, i)] = calc_distances(full_mllkd, mllkds{kk,i}');
        lap0_nbasis(kk, i) = prior_rs{kk,i}.DoF;
    end
end

figure;

subplot(2,1,1)
plot(1:n_iter, lap0_nbasis, 'o-')
xlabel('iteration')
ylabel('basis dimension')
subplot(2,1,2)
semilogy(1:n_iter, lap0_dkl, 'o-')
xlabel('iteration')
ylabel('KL divergence')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load run_iter_zero_ip_bd.mat
pr0_dkl = zeros(n_trial, n_iter);
pr0_dh = zeros(n_trial, n_iter);
pr0_nbasis = zeros(n_trial, n_iter);
for kk = 1:n_trial
    for i = 1:n_iter
        [pr0_dkl(kk, i), pr0_dh(kk, i)] = calc_distances(full_mllkd, mllkds{kk,i}');
        pr0_nbasis(kk, i) = prior_rs{kk, i}.DoF;
    end
end


figure;

subplot(2,1,1)
plot(1:n_iter, pr0_nbasis, 'o-')
xlabel('iteration')
ylabel('basis dimension')
subplot(2,1,2)
semilogy(1:n_iter, pr0_dkl, 'o-')
xlabel('iteration')
ylabel('KL divergence')

%{
load run_iter_lap_ip.mat
lap0_dkl = zeros(n_trial, n_iter);
lap0_dh = zeros(n_trial, n_iter);
lap0_nbasis = zeros(n_trial, n_iter);
for kk = 1:n_trial
    for i = 1:n_iter
        [lap0_dkl(kk, i), lap0_dh(kk, i)] = calc_distances(full_mllkd, mllkds{kk,i}');
        lap0_nbasis(kk, i) = prior_rs{kk,i}.DoF;
    end
end

figure;

subplot(2,1,1)
plot(1:n_iter, lap0_nbasis, 'o-')
xlabel('iteration')
ylabel('basis dimension')
subplot(2,1,2)
semilogy(1:n_iter, lap0_dkl, 'o-')
xlabel('iteration')
ylabel('KL divergence')
%}