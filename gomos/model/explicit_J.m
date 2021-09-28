function Jx = explicit_J(model, HI)

Jx  = -bsxfun(@times, model.AC, HI.d); % jacobian wrt x

% Jx  = -bsxfun(@times,model.AC,exp(-model.AC*HI.x));    

end