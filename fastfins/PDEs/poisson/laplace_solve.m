function soln = laplace_solve(FEM, cond)
% LAPLACE_FOM_SOLVE   solve the linear system of the full FEM by Cholesky
%
% soln = LAPLACE_FOM_SOLVE(FEM, cond)
%
% Tiangang Cui, 02/May/2012

%global forward_count
%forward_count = forward_count + 1;

x           = sparse(1:length(cond),1:length(cond),cond);
A           = FEM.W1*x*FEM.W1' + FEM.W2*x*FEM.W2' + FEM.W3*x*FEM.W3' + FEM.c;

%soln.upper  = chol(A(FEM.p,FEM.p));  % Cholesky factorise 
%soln.lower  = soln.upper'; 

[soln.lower, soln.upper] = lu(A(FEM.p,FEM.p));   % LU

U           = soln.upper\(soln.lower\FEM.fs(FEM.p,:));
soln.G      = full(U(FEM.r,:));
soln.d      = FEM.C*soln.G;

