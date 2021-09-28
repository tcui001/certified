function u = cov_Lv(cov, v)
%COV_LV
%
% v ~ N(0, I)
% u ~ N(0, C)
%
% Tiangang Cui, 17/Jan/2014

u   = cov.RC'*v;

end
