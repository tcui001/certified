function FEM = heat_make_FEM_T1(FEM, mesh)
%HEAT_MAKE_FEM_T1
%
% setup the transient heat problem with Neumann b.c.
%
% Tiangang Cui, 09/May/2014

% tol = 1e-10;

FEM     = heat_make_FEM(FEM, mesh); % setup the basic structure

%{
f1 = 1*exp(-0.5*sum((mesh.node-[0.02;0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2) + ...
     2*exp(-0.5*sum((mesh.node-[max(mesh.gx)-0.02;0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2) + ...
     3*exp(-0.5*sum((mesh.node-[max(mesh.gx)-0.02;max(mesh.gy)-0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2) ;

f2 = 6*exp(-0.5*sum((mesh.node-[0.02;max(mesh.gy)-0.02]*ones(1,mesh.Nnode)).^2)/0.1^2)'/(2*pi*0.1^2);
%}

f1  = exp(-0.5*sum((mesh.node-[0.5;0.5]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
f2  = exp(-0.5*sum((mesh.node-[2.5;0.5]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
f   = f1 - f2;

FEM.fs = (5*1E-5)*(FEM.M*f);
FEM.pe = speye(mesh.Nnode)*0;

end
