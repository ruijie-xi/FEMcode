function mat = local_matrix_FE(i_elem,coe_fun,mesh,mesh_trial,mesh_test,i_vec_trial,dx_trial,dy_trial,i_vec_test,dx_test,dy_test,...
    vec_FE,mesh_FE,i_vec_FE,dx_FE,dy_FE,Gauss_type)
% assemble local matrix

vertices = mesh.P(:,mesh.T(:,i_elem));
[Gauss_weights,Gauss_pts] = generate_Gauss_local_triangle(vertices,Gauss_type);

vec = FE_function_local_2D_read(vec_FE,mesh,mesh_FE,i_elem,i_vec_FE,dx_FE,dy_FE);

basis_test = basis_local_read(mesh_test,i_elem,i_vec_test,dx_test,dy_test);
basis_trial = basis_local_read(mesh_trial,i_elem,i_vec_trial,dx_trial,dy_trial);

mat = vec.*(Gauss_weights.*coe_fun(Gauss_pts).*basis_test)*basis_trial';

