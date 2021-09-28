function Gamma=gomos_prior_GP(alts,sig,L,tau)
% GP prior for gomos example

Gamma1=covfun(alts,sig(1),L(1),tau(1));
Gamma2=covfun(alts,sig(2),L(2),tau(2));
Gamma3=covfun(alts,sig(3),L(3),tau(3));
Gamma4=covfun(alts,sig(4),L(4),tau(4));

Gamma=zeros(length(alts)*4);
Gamma(1:4:end,1:4:end)=Gamma1;
Gamma(2:4:end,2:4:end)=Gamma2;
Gamma(3:4:end,3:4:end)=Gamma3;
Gamma(4:4:end,4:4:end)=Gamma4;

end