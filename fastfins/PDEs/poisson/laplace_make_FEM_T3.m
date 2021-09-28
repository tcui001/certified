function FEM = laplace_make_FEM_T3(FEM, mesh)
%HEAT_MAKE_FEM_T2
%
% Tiangang Cui, 31/Oct/2012

% tol = 1e-10;

FEM = laplace_make_FEM(FEM, mesh);

f1 = 0;
% bottom
for i = 0:12
    f1 = f1 + 1*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(i/12)    ;0                  ]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
end
% top
for i = 0:12
    f1 = f1 + 1*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(i/12)    ; max(mesh.gy)      ]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
end
% count as 12

% left
for i = 1:3
    f1 = f1 + 1*exp(-0.5*sum((mesh.node-[0                      ;max(mesh.gy)*(i/4) ]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
end
% 1
for i = 1:3
    f1 = f1 + 1*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(1/3)     ;max(mesh.gy)*(i/4) ]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
end
% 2
for i = 1:3
    f1 = f1 + 1*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(2/3)     ;max(mesh.gy)*(i/4) ]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
end
% right
for i = 1:3
    f1 = f1 + 1*exp(-0.5*sum((mesh.node-[max(mesh.gx)           ;max(mesh.gy)*(i/4) ]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);
end
% count as 6
 
f2 = 6*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(1/6)  ;max(mesh.gy)/2]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2) + ...
     6*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(3/6)  ;max(mesh.gy)/2]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2) + ...
     6*exp(-0.5*sum((mesh.node-[max(mesh.gx)*(5/6)  ;max(mesh.gy)/2]*ones(1,mesh.Nnode)).^2)/0.05^2)'/(2*pi*0.05^2);

f = f1 - f2;

FEM.fs = FEM.M*f;

end