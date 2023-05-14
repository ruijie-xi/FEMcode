function max_value = Gauss_Linf_error(exact_fun,vector,Gauss_pts,mesh,mesh_FE,i_elem,i_vec)
% compute local Linf error.

fun_value = exact_fun(Gauss_pts,0,0);
fun_value = fun_value(i_vec,:);
value = abs(fun_value-FE_function_local_2D_read(vector,mesh,mesh_FE,i_elem,i_vec,0,0));

max_value = max(value);