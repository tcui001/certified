function Jx = explicit_J(model, HI)
%ADJOINT_JACOBIAN
%
% Tiangang Cui, 19/Mar/2014

switch model.problem
    case{'Laplace'}
        Jx = laplace_jacobian(model, HI);
    case{'Heat'}
        Jx =    heat_jacobian(model, HI);
end

end