function vector = interpolate(fun,mesh,mesh_FE)
% interpolate functions into an FE space.

vector = zeros(mesh_FE.N_node,1);
for i=1:mesh_FE.N_node
    vector(i) = compute_dof(fun,mesh,mesh_FE,i);
end