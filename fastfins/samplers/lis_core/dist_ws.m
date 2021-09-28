function d = dist_ws(param1, param2)
% compute the distance of the Gaussians
% Tiangang Cui, 05/Mar/2014

s1 = (param1.S.^2/sum(param1.S.^2)).^(0.25);
s2 = (param2.S.^2/sum(param2.S.^2)).^(0.25);

tmp = scale_cols(param1.P, s1)'*scale_cols(param2.P,s2);

d = sqrt(1 - sum(tmp(:).^2));

%disp(d);