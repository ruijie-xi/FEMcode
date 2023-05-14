function mat = local_matrix(i_elem,coe_fun,mesh,mesh_trial,mesh_test,i_vec_trial,dx_trial,dy_trial,i_vec_test,dx_test,dy_test,Gauss_type)
% assemble local matrix

vertices = mesh.P(:,mesh.T(:,i_elem));
[Gauss_weights,Gauss_pts] = generate_Gauss_local_triangle(vertices,Gauss_type);

% basis_test = basis_local_2D(Gauss_pts,mesh,mesh_test,i_elem,dx_test,dy_test);
% basis_test = reshape(basis_test(i_vec_test,:,:),mesh_test.N_lb,[]);
basis_test = basis_local_read(mesh_test,i_elem,i_vec_test,dx_test,dy_test);

% basis_trial = basis_local_2D(Gauss_pts,mesh,mesh_trial,i_elem,dx_trial,dy_trial);
% basis_trial = reshape(basis_trial(i_vec_trial,:,:),mesh_trial.N_lb,[]);
basis_trial = basis_local_read(mesh_trial,i_elem,i_vec_trial,dx_trial,dy_trial);

mat = (Gauss_weights.*coe_fun(Gauss_pts).*basis_test)*basis_trial';

