function Jty = matvec_Jty(model, HI, dy)
%MATVEC_JX
%
% Tiangang Cui, 19/Mar/2014

switch model.problem
    case{'Laplace'}
        Jty = laplace_matvec_Jty(model, HI, dy);
    case{'Heat'}
        Jty =    heat_matvec_Jty(model, HI, dy);
    case{'RD'}
        Jty =      RD_matvec_Jty(model, HI, dy);
end

end