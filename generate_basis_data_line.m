function [basis,basis_dx,basis_dy] = generate_basis_data_line(mesh,mesh_FE,Gauss_type,bndy_number)
% compute and cache basis data on edges.

basis = zeros(mesh_FE.dim,3*mesh.N_elem*mesh_FE.N_lb,Gauss_type);
basis_dx = zeros(mesh_FE.dim,3*mesh.N_elem*mesh_FE.N_lb,Gauss_type);
basis_dy = zeros(mesh_FE.dim,3*mesh.N_elem*mesh_FE.N_lb,Gauss_type);

if isempty(bndy_number)
    flag = true(1,mesh.N_edge);
else
    flag = false(1,mesh.N_edge);
    for i = 1:mesh.N_edge
        if mesh.E(1,i)==bndy_number
            flag(i) = true;
        end
    end
end

for i_elem = 1:mesh.N_elem
    for i = 1:mesh.type
        i_edge = abs(mesh.TE(i,i_elem));
        if flag(i_edge)
            endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
    
            [~,pts] = generate_Gauss_local_line(endpts,Gauss_type);
            range = (i_elem-1)*3*mesh_FE.N_lb+(i-1)*mesh_FE.N_lb+(1:mesh_FE.N_lb);
            basis(:,range,:) = basis_local_2D(pts,mesh,mesh_FE,i_elem,0,0);
            basis_dx(:,range,:) = basis_local_2D(pts,mesh,mesh_FE,i_elem,1,0);
            basis_dy(:,range,:) = basis_local_2D(pts,mesh,mesh_FE,i_elem,0,1);
        end
    end
end