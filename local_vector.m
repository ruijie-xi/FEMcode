function vec = local_vector(i_elem,coe_fun,mesh,mesh_test,...
    i_vec_test,dx_test,dy_test,Gauss_type)
% assemble local vector

vertices = mesh.P(:,mesh.T(:,i_elem));
[Gauss_weights,Gauss_pts] = generate_Gauss_local_triangle(vertices,Gauss_type);

% basis_test = basis_local_2D(Gauss_pts,mesh,mesh_test,i_elem,dx_test,dy_test);
% basis_test = reshape(basis_test(i_vec_test,:,:),mesh_test.N_lb,[]);
basis_test = basis_local_read(mesh_test,i_elem,i_vec_test,dx_test,dy_test);

vec = basis_test*(Gauss_weights.*coe_fun(Gauss_pts))';


