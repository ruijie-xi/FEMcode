function plot_boundary(mesh,bndy_number)
% plot boundary pts.

x = [];
y = [];
for i_edge = 1:mesh.N_edge
    if mesh.E(1,i_edge)==bndy_number
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        x = [x,endpts(1,:)];
        y = [y,endpts(2,:)];
    end
end
plot(x,y,'ro')