function w  = matvec_dGds(FEM, s, ds, dt)
%MATVEC_DGDS
%
% Tiangang Cui, 20/May/2014

if FEM.ML_flag
    w   = FEM.M*ds + dt*(FEM.K*ds)    - dt*FEM.a*(FEM.M*ds) + 2*dt*FEM.a*(FEM.M*(s.*ds));
else    
    w   = ds       + dt*(FEM.iMLK*ds) - dt*FEM.a*ds         + 2*dt*FEM.a*(s.*ds);
end

end