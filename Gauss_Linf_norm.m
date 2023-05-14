function max_value = Gauss_Linf_norm(vector,Gauss_pts,mesh,mesh_FE,i_elem,i_vec)
% compute local Linf norm.

value = abs(FE_function_local_2D_read(vector,mesh,mesh_FE,i_elem,i_vec,0,0));

max_value = max(value);