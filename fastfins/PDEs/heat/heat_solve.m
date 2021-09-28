function soln = heat_solve(FEM, x)
%HEAT_SOLVE  
%
% solves the transient heat equation by Cholesky and implicit Euler
%
% Tiangang Cui, 09/May/2014

% plot(cond)

setup.type      = 'ict';
setup.michol    = 'on';
setup.droptol   = 0.1;
setup.shape     = 'lower';

TOL         = 1E-6;
MAXIT       = 50;

K           = sparse(1:length(x),1:length(x),x);
soln.N      = FEM.W1*K*FEM.W1' + FEM.W2*K*FEM.W2' + FEM.W3*K*FEM.W3' + FEM.pe;

soln.G      = zeros(FEM.DoF, FEM.Nsteps+1);
soln.L      = cell(FEM.Nsteps+1,1);

for i = 1:FEM.Nsteps
    g       = FEM.M*soln.G(:,i) + FEM.Tsteps(i)*FEM.fs; % fs is the integrated forcing over time FEM.ts(i)
    B       = FEM.M + FEM.Tsteps(i)*soln.N;
    
    if FEM.GMRES_flag
        soln.L{i+1}         = ichol(B, setup);
        [soln.G(:,i+1),~]   = minres(B, g, TOL, MAXIT, soln.L{i+1}, soln.L{i+1}');
    else
        soln.L{i+1}         = chol(B(FEM.p,FEM.p))';
        soln.G(FEM.p,i+1)   = soln.L{i+1}'\(soln.L{i+1}\g(FEM.p));
    end
end

soln.d      = (FEM.C*soln.G)*FEM.obs_int;
soln.NoP    = size(x, 1);

end

%{
size(soln.d)

d = zeros(size(soln.d));
for i = 1:FEM.max_Nt
    d(:,FEM.obs_ind{i}) = soln.d_pri(:,i)*FEM.obs_w0{i} + soln.d_pri(:,i+1)*FEM.obs_w1{i};
end

sum(d(:)-soln.d(:))
%}