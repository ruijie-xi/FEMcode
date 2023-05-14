function int_value = Gauss_quad_2D_norm(vector,Gauss_pts,Gauss_weights,...
    mesh,mesh_FE,i_elem,dx,dy,i_vec)
% compute local norm.

int_value = (FE_function_local_2D_read(vector,mesh,mesh_FE,i_elem,i_vec,dx,dy)).^2*Gauss_weights';