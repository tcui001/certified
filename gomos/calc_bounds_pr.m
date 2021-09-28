function bounds = calc_bounds_pr(ref_H, prior, nps)

[V, ~] = cov_eig(prior.cov);

U1 = matvec_prior_Lt(prior, V);
U2 = matvec_prior_invL(prior, V);

m = length(nps);
bounds = zeros(m, 1);

tH = trace(ref_H);

for i = 1:m
    ind = 1:nps(i);
    a = trace(U1(:,ind)'*(ref_H*U2(:,ind)));
    b = trace( (U2(:,ind)'*(ref_H*U2(:,ind))) * (U1(:,ind)'*U1(:,ind)) );
    bounds(i) = tH - 2*a + b;
end

end