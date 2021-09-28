function Jx = matvec_Jx(model, HI, dx)
%MATVEC_JX
%
% Tiangang Cui, 19/Mar/2014

switch model.problem
    case{'Laplace'}
        Jx  = laplace_matvec_Jx(model, HI, dx);
    case{'Heat'}
        Jx  =    heat_matvec_Jx(model, HI, dx);
    case{'RD'}
        Jx  =      RD_matvec_Jx(model, HI, dx);
end

end