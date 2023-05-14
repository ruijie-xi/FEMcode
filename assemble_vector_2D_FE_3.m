function b = assemble_vector_2D_FE_3(rhs_fun,coe_num,mesh,mesh_test,...
    vec_FE_1,mesh_FE_1,vec_FE_2,mesh_FE_2,vec_FE_3,mesh_FE_3,vec_info,Gauss_type)
% assemble vector with 3 FE functions.

vector_size = mesh_test.N_node;
num_vec = size(vec_info,1);
local_vec_size = mesh_test.N_lb;
index_i = zeros(mesh.N_elem*mesh_test.N_lb*num_vec,1);
index_j = ones(mesh.N_elem*mesh_test.N_lb*num_vec,1);
index_v = zeros(mesh.N_elem*mesh_test.N_lb*num_vec,1);

for i_vec = 1:num_vec
    i_vec_test = vec_info(i_vec,1);
    dx_test = vec_info(i_vec,2);
    dy_test = vec_info(i_vec,3);
    i_vec_FE_1 = vec_info(i_vec,4);
    dx_FE_1 = vec_info(i_vec,5);
    dy_FE_1 = vec_info(i_vec,6);
    i_vec_FE_2 = vec_info(i_vec,7);
    dx_FE_2 = vec_info(i_vec,8);
    dy_FE_2 = vec_info(i_vec,9);
    i_vec_FE_3 = vec_info(i_vec,10);
    dx_FE_3 = vec_info(i_vec,11);
    dy_FE_3 = vec_info(i_vec,12);

    for n = 1:mesh.N_elem
        vec = coe_num(i_vec)*local_vector_FE_3(n,rhs_fun,mesh,mesh_test,i_vec_test,dx_test,dy_test,...
            vec_FE_1,mesh_FE_1,i_vec_FE_1,dx_FE_1,dy_FE_1,vec_FE_2,mesh_FE_2,i_vec_FE_2,dx_FE_2,dy_FE_2,...
            vec_FE_3,mesh_FE_3,i_vec_FE_3,dx_FE_3,dy_FE_3,Gauss_type);
        i = mesh_test.T((i_vec_test-1)*mesh_test.N_lb+1:i_vec_test*mesh_test.N_lb,n);
        range = (i_vec-1)*local_vec_size*mesh.N_elem+(n-1)*mesh_test.N_lb+(1:mesh_test.N_lb);
        index_i(range) = i;
        index_v(range) = vec;
    end
end

b = sparse(index_i,index_j,index_v,vector_size,1);