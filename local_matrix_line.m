function mat = local_matrix_line(i_edge,coe_fun,mesh,mesh_trial,mesh_test,i_vec_trial,dx_trial,dy_trial,i_vec_test,dx_test,dy_test,Gauss_type)
% assemble local matrix (integral on line)

i_elem = mesh.E(2,i_edge);
endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
[Gauss_weights,Gauss_pts] = generate_Gauss_local_line(endpts,Gauss_type);

mat1 = Gauss_weights.*coe_fun(Gauss_pts).*...
    basis_local_read_line(mesh,mesh_test,i_elem,i_edge,i_vec_test,dx_test,dy_test);

mat = mat1*basis_local_read_line(mesh,mesh_trial,i_elem,i_edge,i_vec_trial,dx_trial,dy_trial)';


