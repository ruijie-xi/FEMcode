function result = FE_function_local_2D(coords,vector,mesh,mesh_FE,i_elem,i_vec,dx,dy)
% evaluate an FE function at any pt in a certain element.

local_vec = get_local_dof(vector,mesh,mesh_FE,i_elem,i_vec);

basis = basis_local_2D(coords,mesh,mesh_FE,i_elem,dx,dy);
basis = reshape(basis(i_vec,:,:),mesh_FE.N_lb,[]);

result = local_vec'*basis;
