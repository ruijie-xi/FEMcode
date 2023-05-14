function mat = local_matrix_FE_FE(i_elem,coe_fun,mesh,mesh_trial,mesh_test,i_vec_trial,dx_trial,dy_trial,i_vec_test,dx_test,dy_test,...
    vec_FE_1,mesh_FE_1,i_vec_FE_1,dx_FE_1,dy_FE_1,vec_FE_2,mesh_FE_2,i_vec_FE_2,dx_FE_2,dy_FE_2,Gauss_type)
% assemble local matrix

vertices = mesh.P(:,mesh.T(:,i_elem));
[Gauss_weights,Gauss_pts] = generate_Gauss_local_triangle(vertices,Gauss_type);

vec_1 = FE_function_local_2D_read(vec_FE_1,mesh,mesh_FE_1,i_elem,i_vec_FE_1,dx_FE_1,dy_FE_1);
vec_2 = FE_function_local_2D_read(vec_FE_2,mesh,mesh_FE_2,i_elem,i_vec_FE_2,dx_FE_2,dy_FE_2);


basis_test = basis_local_read(mesh_test,i_elem,i_vec_test,dx_test,dy_test);
basis_trial = basis_local_read(mesh_trial,i_elem,i_vec_trial,dx_trial,dy_trial);

mat = vec_1.*vec_2.*(Gauss_weights.*coe_fun(Gauss_pts).*basis_test)*basis_trial';

