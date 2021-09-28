function FEM = laplace_make_FEM_DILI(FEM, mesh)
%HEAT_MAKE_FEM_T1
%
% Tiangang Cui, 31/Oct/2012

% tol = 1e-10;

FEM = laplace_make_FEM(FEM, mesh);

% put Gaussian source in the middle, and constant flux around the boundary
% Gaussian force
f = 2*exp(-0.5*sum((mesh.node-[0.3;0.3]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2) + ...
    3*exp(-0.5*sum((mesh.node-[0.7;0.7]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2) - ...
    2*exp(-0.5*sum((mesh.node-[0.3;0.7]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2) - ...
    3*exp(-0.5*sum((mesh.node-[0.7;0.3]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);

FEM.fs = FEM.M*f;

end