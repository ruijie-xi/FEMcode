function b = assemble_vector_2D(rhs_fun,mesh,mesh_test,dx_test,dy_test,i_vec_test,Gauss_type)
% assemble vector like \int f ((\partial^m_x\partial^n_y v_l)) dxdy.
% where:
% - f = rhs_fun
% - m = dx_test
% - n = dy_test
% - l = i_vec_test

vector_size = mesh_test.N_node;
index_i = zeros(mesh.N_elem*mesh_test.N_lb,1);
index_j = ones(mesh.N_elem*mesh_test.N_lb,1);
index_v = zeros(mesh.N_elem*mesh_test.N_lb,1);

for n = 1:mesh.N_elem
    vec = local_vector(n,rhs_fun,mesh,mesh_test,i_vec_test,dx_test,dy_test,Gauss_type);
    i = mesh_test.T((i_vec_test-1)*mesh_test.N_lb+1:i_vec_test*mesh_test.N_lb,n);
    range = (n-1)*mesh_test.N_lb+(1:mesh_test.N_lb);
    index_i(range) = i;
    index_v(range) = vec;
end

b = sparse(index_i,index_j,index_v,vector_size,1);