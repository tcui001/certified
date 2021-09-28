function plot_marginal(samples, ri, c, subplot)

if length(ri) > 3
    ii = ri(3:end);
    n = length(ii);
else
    n = ri(3);
    ii = 1:n;
end

for i = 1:n
    j = ii(i);
    subplot(ri(1), ri(2), i);
    [f, xi] = ksdensity(samples(j,:));
    plot(xi, f, c, 'linewidth', 2, 'markersize', 1); hold on
    % title(['KL mode ' num2str(j)])
    xlabel(['$x_' num2str(j) '$'])
end


end
