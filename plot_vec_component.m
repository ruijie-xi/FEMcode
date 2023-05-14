function plot_vec_component(vector,mesh,mesh_FE,component)

X = zeros(mesh.type,mesh.N_elem);
Y = zeros(mesh.type,mesh.N_elem);
U = zeros(mesh.type,mesh.N_elem);
V = zeros(mesh.type,mesh.N_elem);

for i_elem = 1:mesh.N_elem
    X(:,i_elem) = mesh.P(1,mesh.T(:,i_elem));
    Y(:,i_elem) = mesh.P(2,mesh.T(:,i_elem));
    
    vertices = mesh.P(:,mesh.T(:,i_elem));
    if component==1
        U(:,i_elem) = FE_function_local_2D(vertices,vector,mesh,mesh_FE,i_elem,1,0,0);
    elseif component==2
        V(:,i_elem) = FE_function_local_2D(vertices,vector,mesh,mesh_FE,i_elem,2,0,0);
    end
end

if component==1
    patch(X,Y,U,'LineStyle','none');
    colorbar
elseif component==2
    patch(X,Y,V,'LineStyle','none');
    colorbar

end
