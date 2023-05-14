function int_value = Gauss_quad_2D_error(exact_fun,vector,Gauss_pts,Gauss_weights,mesh,mesh_FE,i_elem,dx,dy,i_vec)
% compute local L2 error.

fun_value = exact_fun(Gauss_pts,dx,dy);
fun_value = fun_value(i_vec,:);
int_value = (fun_value-FE_function_local_2D_read(vector,mesh,mesh_FE,i_elem,i_vec,dx,dy)).^2*Gauss_weights';