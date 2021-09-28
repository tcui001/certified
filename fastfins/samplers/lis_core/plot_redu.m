function plot_redu(prior, param, m, n)
% plot the basis of a given param struct
% Tiangang Cui, 01/May/2013

tmp = prior_L_mult(prior, param.P);
j = 30;
for i = 1:size(tmp,2)
    if j == 30
        figure('position',[100 100 750 1000])
        j = 0;
    end
    j = j+1;
    subplot(6,5,j)
    imagesc(reshape(tmp(:,i),m,n));
    title(['Modes ' num2str(i)])
end
