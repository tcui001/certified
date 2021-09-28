function v = cov_invLu(cov, u)
%COV_INVLU 
%
% whitening transformation
%
% Tiangang Cui, 17/Jan/2014

v   = cov.RC'\u;

end