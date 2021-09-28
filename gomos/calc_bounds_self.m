function bounds = calc_bounds_self(H, nps)

[~, D] = svd(H);
d = diag(D);

m = length(nps);
bounds = zeros(m, 1);


for i = 1:m
    ind = 1:nps(i);
    bounds(i) = sum(d) - sum(d(ind));
end

end