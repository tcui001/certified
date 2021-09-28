
function [dkl, dh] = calc_distances(full_mllkd, in_mllkd)

n = length(full_mllkd);
nrp = size(in_mllkd, 1);

f_mllkd = repmat(full_mllkd, nrp, 1);

t = f_mllkd - in_mllkd;

%{
logratio2 = zeros(nrp, 1);
for i = 1:nrp
    mt = max(t(i,:));
    logratio2(i) = mt + log(sum(exp(t(i,:)-mt), 2)) - log(n);
end
%}

logratio = logsumexp(t);

%logratio - logratio2

% dkl
dkl = ( sum(-t, 2)/n + logratio )';

% dh

dh = sqrt(1 - exp( logsumexp(0.5*t) - 0.5*logratio ))';

%dh = 1 - exp(- logratio) .*  ( sum(exp(0.5 * t), 2) /n );

end

function a = logsumexp(t)

[m, n] = size(t);
a = zeros(m, 1);

for i = 1:m
    mt   = max(t(i,:));
    a(i) = mt + log(sum(exp(t(i,:)-mt), 2)) - log(n);
end

end