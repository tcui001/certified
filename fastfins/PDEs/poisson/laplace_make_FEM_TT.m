function FEM = laplace_make_FEM_TT(FEM, mesh)
%HEAT_MAKE_FEM_T1
%
% Tiangang Cui, 31/Oct/2012

% tol = 1e-10;

tol = 0.05;

FEM = laplace_make_FEM(FEM, mesh);

f1 = 2*exp(-0.5*sum((mesh.node-[max(mesh.gx)-0.02;0.02]*ones(1,mesh.Nnode)).^2)/tol^2)'/(2*pi*tol^2) + ...
     4*exp(-0.5*sum((mesh.node-[max(mesh.gx)-0.02;max(mesh.gy)-0.02]*ones(1,mesh.Nnode)).^2)/tol^2)'/(2*pi*tol^2) ;

f2 = 3*exp(-0.5*sum((mesh.node-[0.02;0.02]*ones(1,mesh.Nnode)).^2)/tol^2)'/(2*pi*tol^2) + ...
     3*exp(-0.5*sum((mesh.node-[0.02;max(mesh.gy)-0.02]*ones(1,mesh.Nnode)).^2)/tol^2)'/(2*pi*tol^2);

f = (f1 - f2/sum(f2)*sum(f1));

FEM.fs = 1E-5*(FEM.M*f);

FEM.c = FEM.Mb*(mesh.Nside*mesh.Mside);

end
