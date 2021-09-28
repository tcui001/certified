function FEM = heat_make_FEM_T2(FEM, mesh)
%HEAT_MAKE_FEM_T1
%
% setup the transient heat problem with Neumann b.c.
%
% Tiangang Cui, 09/May/2014

tol = 1e-10;

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

FEM.f  = (5*1E-5)*(FEM.M*f);

% penalty term for the boundary
pe_t    = sparse([],[],[], mesh.Nnode, mesh.Nnode, 4*mesh.Nbndf);
pe_b    = sparse([],[],[], mesh.Nnode, mesh.Nnode, 4*mesh.Nbndf);
fb      = sparse([],[],[], mesh.Nnode, 1, mesh.Nbndf);
ft      = sparse([],[],[], mesh.Nnode, 1, mesh.Nbndf);

% boundary
yb      = min(mesh.node(2,:));
yt      = max(mesh.node(2,:));


for i = 1:mesh.Nbndf
    ind                 = mesh.node_map_bnd(:,i);
    dx                  = norm( mesh.node(:,ind(2)) - mesh.node(:,ind(1)) );
    
    % bottom boundary 
    if abs( mean(mesh.node(2,ind)) - yb ) < tol 
        pe_b(ind,ind)   = pe_b(ind,ind) + dx*mesh.locmass_bnd;
        fb(ind)         = 1;
    end
    
    % top boundary
    if abs( mean(mesh.node(2,ind)) - yt ) < tol 
        pe_t(ind,ind)   = pe_t(ind,ind) + dx*mesh.locmass_bnd;
        ft(ind)         = 1;
    end
end

FEM.pe  = (pe_b + pe_t)*(mesh.Nside*mesh.Mside);
FEM.ft  = FEM.pe*ft;
FEM.fb  = FEM.pe*fb;

end
