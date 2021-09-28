function movie_heat(mesh, FEM, soln)
% plot
figure('position',[100 100 700 300])

temp = cumsum(FEM.ts);

for i = 2:(FEM.max_Nt+1)
    pause(0.5)
    subplot(1,2,1)
    meshc(mesh.gx,mesh.gx,reshape(soln.G(:,i),mesh.N_side+1,mesh.N_side+1));
    title(['Solution at time ' num2str(temp(i-1))])
    xlabel('x')
    ylabel('y')
    
    subplot(1,2,2)
    meshc(mesh.gx,mesh.gx,reshape(soln.G(:,i)-soln.G(:,i-1),mesh.N_side+1,mesh.N_side+1));
    title('Difference of Solutions')
    xlabel('x')
    ylabel('y')
end
