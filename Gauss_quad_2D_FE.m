function int_value = Gauss_quad_2D_FE(vector,Gauss_pts,Gauss_weights,...
    mesh,mesh_FE,i_elem,dx,dy,i_vec)
% compute local integral value.

int_value = FE_function_local_2D_read(vector,mesh,mesh_FE,i_elem,i_vec,dx,dy)*Gauss_weights';