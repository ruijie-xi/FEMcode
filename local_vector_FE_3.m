function vec = local_vector_FE_3(i_elem,coe_fun,mesh,mesh_test,...
    i_vec_test,dx_test,dy_test,vec_FE_1,mesh_FE_1,i_vec_FE_1,dx_FE_1,dy_FE_1,...
    vec_FE_2,mesh_FE_2,i_vec_FE_2,dx_FE_2,dy_FE_2,...
    vec_FE_3,mesh_FE_3,i_vec_FE_3,dx_FE_3,dy_FE_3,Gauss_type)
% assemble local vector

vertices = mesh.P(:,mesh.T(:,i_elem));
[Gauss_weights,Gauss_pts] = generate_Gauss_local_triangle(vertices,Gauss_type);

vec_1 = FE_function_local_2D_read(vec_FE_1,mesh,mesh_FE_1,i_elem,i_vec_FE_1,dx_FE_1,dy_FE_1);
vec_2 = FE_function_local_2D_read(vec_FE_2,mesh,mesh_FE_2,i_elem,i_vec_FE_2,dx_FE_2,dy_FE_2);
vec_3 = FE_function_local_2D_read(vec_FE_3,mesh,mesh_FE_3,i_elem,i_vec_FE_3,dx_FE_3,dy_FE_3);

% basis_test = basis_local_2D(Gauss_pts,mesh,mesh_test,i_elem,dx_test,dy_test);
% basis_test = reshape(basis_test(i_vec_test,:,:),mesh_test.N_lb,[]);
basis_test = basis_local_read(mesh_test,i_elem,i_vec_test,dx_test,dy_test);

vec = vec_1.*vec_2.*vec_3.*basis_test*(Gauss_weights.*coe_fun(Gauss_pts))';