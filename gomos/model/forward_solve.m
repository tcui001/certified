function soln = forward_solve(model, x)

d       = exp(-model.C*reshape(x,model.ngas,model.nalts)*model.A'); % forward solve
soln.d  = d(:);
%soln.x  = x;

end