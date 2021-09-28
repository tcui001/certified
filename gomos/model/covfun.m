function C=covfun(x,th,L,tau)

n=length(x);
C=zeros(n);

for i=1:n
    for j=i:n
        C(i,j)=th*exp(-(x(i)-x(j))^2/2/L);
        C(j,i)=C(i,j);
        if i==j
            C(i,j)=C(i,j)+tau;
        end
    end
end