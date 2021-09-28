function soln = forward_solve(model, x)
%FORWARD_SOLVE  
%
% A driver function of the forward solver
%
% Tiangang Cui, 10/Mar/2013

switch model.problem
    case{'Laplace'}
        soln = laplace_solve(model, x);
    case{'Heat'}
        soln =    heat_solve(model, x);
    case{'RD'}
        soln =      RD_solve(model, x);
end

end