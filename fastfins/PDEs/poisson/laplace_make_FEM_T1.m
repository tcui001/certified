function FEM = laplace_make_FEM_T1(FEM, mesh)
%HEAT_MAKE_FEM_T1
%
% Tiangang Cui, 31/Oct/2012

% tol = 1e-10;

FEM = laplace_make_FEM(FEM, mesh);

f1 = 1*exp(-0.5*sum((mesh.node-[0.02;0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2) + ...
     2*exp(-0.5*sum((mesh.node-[max(mesh.gx)-0.02;0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2) + ...
     3*exp(-0.5*sum((mesh.node-[max(mesh.gx)-0.02;max(mesh.gy)-0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2) ;

f2 = 6*exp(-0.5*sum((mesh.node-[0.02;max(mesh.gy)-0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2);

f = (f1 - f2/sum(f2)*sum(f1));

FEM.fs = FEM.M*f;

FEM.c = FEM.Mb*(mesh.Nside*mesh.Mside);

end
