function plot_scalar(vector,mesh,mesh_FE)
% plot a scalar function.

clf;
X = zeros(mesh.type,mesh.N_elem);
Y = zeros(mesh.type,mesh.N_elem);
Z = zeros(mesh.type,mesh.N_elem);

for i_elem = 1:mesh.N_elem
    X(:,i_elem) = mesh.P(1,mesh.T(:,i_elem));
    Y(:,i_elem) = mesh.P(2,mesh.T(:,i_elem));
    
    vertices = mesh.P(:,mesh.T(:,i_elem));
    Z(:,i_elem) = FE_function_local_2D(vertices,vector,mesh,mesh_FE,i_elem,1,0,0);
end
patch(X,Y,Z,'LineStyle','none');
colorbar