function Jx = matvec_Jx(model, HI, dx)    

tmp = - model.C * reshape(dx, model.ngas, model.nalts) * model.A';
Jx  =   tmp(:).*HI.d;

end