function bounds = calc_bounds(ref_H, in_H, nps)

[V, ~] = svd(in_H);

m = length(nps);
bounds = zeros(m, 1);

tH = trace(ref_H);

for i = 1:m
    ind = 1:nps(i);
    bounds(i) = tH - trace( V(:,ind)'*(ref_H*V(:,ind)) );
end

end