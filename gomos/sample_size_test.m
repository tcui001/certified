
n_samples = [50, 100, 200, 500, 1000, 2000];

nps = 1:50;
nrp = length(nps);

self_bound = calc_bounds_self(new_H, nps);

true_bounds = zeros(length(n_samples),nrp);
int_bounds  = zeros(length(n_samples),nrp);

ind = randperm(out_dili.size, n_samples(end)) ;
for kk = 1:length(n_samples)
    
    Hl_g = zeros(prior.DoF);
    for j = ind(1:n_samples(kk))
        [~,~,~,gmllkd] = minus_log_post(model, obs, prior, out_dili.v_samples(j,:)');
        Hl_g = Hl_g + gmllkd*gmllkd';
    end
    
    Hl_g = Hl_g/n_samples(kk);
    
    true_bounds(kk,:)   = calc_bounds(new_H, Hl_g, nps);
    int_bounds(kk,:)    = calc_bounds_self(Hl_g, nps);
    
end

figure('position',  [100 100 900 350])
subplot(1,2,1)
semilogy(nps, self_bound, 'k-', 'linewidth', 2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
semilogy(nps, true_bounds)
ax.ColorOrderIndex = 1;
semilogy(nps, int_bounds, '--')
semilogy(nps, self_bound, 'k-', 'linewidth', 2)
xlabel('rank of the projector')
title('single trial')
h = legend('K=10^6', 'K=50', 'K=100', 'K=200', 'K=500', 'K=1000', 'K=2000', 'location', 'SouthWest' );
set(h, 'box', 'off');
axis([0, 50, 1E-10 1E7])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_samples = [50, 200, 1000];

true_bounds = zeros(length(n_samples),nrp);
int_bounds  = zeros(length(n_samples),nrp);
subplot(1,2,2)
semilogy(nps, self_bound, 'k-', 'linewidth', 2)
hold on
ax = gca;
for ii = 1:10
    
    ind = randperm(out_dili.size, n_samples(end)) ;
    for kk = 1:length(n_samples)
        
        Hl_g = zeros(prior.DoF);
        for j = ind(1:n_samples(kk))
            [~,~,~,gmllkd] = minus_log_post(model, obs, prior, out_dili.v_samples(j,:)');
            Hl_g = Hl_g + gmllkd*gmllkd';
        end
        
        Hl_g = Hl_g/n_samples(kk);
        
        true_bounds(kk,:)   = calc_bounds(new_H, Hl_g, nps);
        int_bounds(kk,:)    = calc_bounds_self(Hl_g, nps);
        
    end
    
    ax.ColorOrderIndex = 1;
    semilogy(nps, true_bounds(1,:))
    ax.ColorOrderIndex = 2;
    semilogy(nps, true_bounds(2,:))
    ax.ColorOrderIndex = 3;
    semilogy(nps, true_bounds(3,:))
    if ii == 1
        h = legend('K=10^6', 'K=50',' K=200', 'K=1000', 'location', 'SouthWest');
    end
    
    ax.ColorOrderIndex = 1;
    semilogy(nps, int_bounds(1,:), '--')
    ax.ColorOrderIndex = 2;
    semilogy(nps, int_bounds(2,:), '--')
    
end

semilogy(nps, self_bound, 'k-', 'linewidth', 2)
h.String = {'K=10^6', 'K=50',' K=200', 'K=1000'};
set(h, 'box', 'off');
xlabel('rank of the projector')
title('multiple trials')
axis([0, 50, 1E-10 1E7])
