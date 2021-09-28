function autocorr_comp(full_samples, outEr, outE2, maxn, lags, indices, ri)


for i = 1:ri(3)
    subplot(ri(1), ri(2), i);
    j = indices(i);
    [af,lf] = autocorr(full_samples(j, 1E5:end), maxn/lags(1));
    [ar,lr] = autocorr(outEr(j, 1E5:end), maxn/lags(2));
    [a2,l2] = autocorr(outE2(j, 1E5:end), maxn/lags(3));
    plot(lf*lags(1), af, 'k', lr*lags(2), ar, 'b:', l2*lags(3), a2, 'r-.')
    title(['KL mode ' num2str(j)])
    set(gca, 'ylim', [0 1])
end

end