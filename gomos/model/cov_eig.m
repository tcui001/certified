function [V, d] = cov_eig(cov)
%COV_EIG
%
% Tiangang Cui, 17/Jan/2014

[V, D]      = eig(cov.C);
[d, ind]    = sort( diag(D), 'descend' );
 V          = V(:, ind);

end
