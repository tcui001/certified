function Jty = matvec_Jty(model, HI, dy)    

Y   = reshape(HI.d(:).*dy(:), size(model.C,1), size(model.A,1));
Jty = - model.C' * Y * model.A;
Jty = Jty(:);

end