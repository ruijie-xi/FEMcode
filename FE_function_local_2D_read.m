function result = FE_function_local_2D_read(vector,mesh,mesh_FE,i_elem,i_vec,dx,dy)
% evaluate an FE function at cached gaussian pts in a certain element.


local_vec = get_local_dof(vector,mesh,mesh_FE,i_elem,i_vec);

basis = basis_local_read(mesh_FE,i_elem,i_vec,dx,dy);

result = local_vec'*basis;
