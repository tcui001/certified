function v = cov_invLtu(cov, u)
%COV_INVLTU
%
% Tiangang Cui, 17/Jan/2014

v   = cov.RC\u;

end