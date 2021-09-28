function soln = heat_solve_simple(FEM, x)
%HEAT_SOLVE_SIMPLE  
%
% solves the transient heat equation by Cholesky and implicit Euler
%
% Tiangang Cui, 09/May/2014

% plot(cond)

K           = sparse(1:length(x),1:length(x),x);
N           = FEM.W1*K*FEM.W1' + FEM.W2*K*FEM.W2' + FEM.W3*K*FEM.W3';

soln.G      = zeros(FEM.DoF, FEM.Nsteps+1);

for i       = 1:FEM.Nsteps
    %
    g       = M*soln.G(:,i) + FEM.fs(i); % fs is the integrated forcing over time FEM.ts(i)
    B       = M + FEM.Tsteps(i)*N;
    R       = chol(B(FEM.p,FEM.p));
    soln.G  (FEM.p,i+1) = R\(R'\g(FEM.p));
end

soln.d_pri  = FEM.C*soln.G;
soln.d      = soln.d_pri*FEM.obs_int;

end

%{
size(soln.d)

d = zeros(size(soln.d));
for i = 1:FEM.max_Nt
    d(:,FEM.obs_ind{i}) = soln.d_pri(:,i)*FEM.obs_w0{i} + soln.d_pri(:,i+1)*FEM.obs_w1{i};
end

sum(d(:)-soln.d(:))
%}