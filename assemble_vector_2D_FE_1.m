function b = assemble_vector_2D_FE_1(rhs_fun,coe_num,mesh,mesh_test,...
    vec_FE,mesh_FE,vec_info,Gauss_type)
% assemble vector with an FE function.
% like like c \int f ((\partial^m_x\partial^n_y v_l)) ((\partial^a_x\partial^b_y w_p) dxdy.
% Each row of vec_info is [l,m,n,p,a,b]

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
    i_vec_FE = vec_info(i_vec,4);
    dx_FE = vec_info(i_vec,5);
    dy_FE = vec_info(i_vec,6);
    for n = 1:mesh.N_elem
        vec = coe_num(i_vec)*local_vector_FE(n,rhs_fun,mesh,mesh_test,i_vec_test,dx_test,dy_test,...
            vec_FE,mesh_FE,i_vec_FE,dx_FE,dy_FE,Gauss_type);
        i = mesh_test.T((i_vec_test-1)*mesh_test.N_lb+1:i_vec_test*mesh_test.N_lb,n);
        range = (i_vec-1)*local_vec_size*mesh.N_elem+(n-1)*mesh_test.N_lb+(1:mesh_test.N_lb);
        index_i(range) = i;
        index_v(range) = vec;
    end
end

b = sparse(index_i,index_j,index_v,vector_size,1);