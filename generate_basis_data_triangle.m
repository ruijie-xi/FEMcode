function [basis,basis_dx,basis_dy] = generate_basis_data_triangle(mesh,mesh_FE,Gauss_type)
% compute and cache basis data on elements.

basis = zeros(mesh_FE.dim,mesh.N_elem*mesh_FE.N_lb,Gauss_type);
basis_dx = zeros(mesh_FE.dim,mesh.N_elem*mesh_FE.N_lb,Gauss_type);
basis_dy = zeros(mesh_FE.dim,mesh.N_elem*mesh_FE.N_lb,Gauss_type);

for i_elem = 1:mesh.N_elem
    vertices = mesh.P(:,mesh.T(:,i_elem));
    [~,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
    range = (i_elem-1)*mesh_FE.N_lb+(1:mesh_FE.N_lb);
    basis(:,range,:) = basis_local_2D(pts,mesh,mesh_FE,i_elem,0,0);
    basis_dx(:,range,:) = basis_local_2D(pts,mesh,mesh_FE,i_elem,1,0);
    basis_dy(:,range,:) = basis_local_2D(pts,mesh,mesh_FE,i_elem,0,1);
end