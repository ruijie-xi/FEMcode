function A = assemble_matrix_2D(coe_fun,coe_num,mesh,mesh_trial,mesh_test,...
    mat_info,Gauss_type)
% assemble matrix like c \int f (\partial^i_x\partial^j_y u_k)
% ((\partial^m_x\partial^n_y v_l)) dxdy.
% where:
% - f = coe_fun
% - c = coe_num
% - each row in mat_info is [k,l,i,j,m,n]

matrix_size = [mesh_test.N_node,mesh_trial.N_node];
num_mat = size(mat_info,1);

local_mat_size = mesh_trial.N_lb*mesh_test.N_lb;
index_i = zeros(local_mat_size*mesh.N_elem*num_mat,1);
index_j = zeros(local_mat_size*mesh.N_elem*num_mat,1);
index_v = zeros(local_mat_size*mesh.N_elem*num_mat,1);

for n = 1:mesh.N_elem

    for i_mat = 1:num_mat

        i_vec_trial = mat_info(i_mat,1);
        i_vec_test = mat_info(i_mat,2);
        dx_trial = mat_info(i_mat,3);
        dy_trial = mat_info(i_mat,4);
        dx_test = mat_info(i_mat,5);
        dy_test = mat_info(i_mat,6);
        
  
        mat = coe_num(i_mat)*local_matrix(n,coe_fun,mesh,mesh_trial,mesh_test,i_vec_trial,dx_trial,dy_trial,i_vec_test,dx_test,dy_test,Gauss_type);
    
        count = (i_mat-1)*local_mat_size*mesh.N_elem+(n-1)*local_mat_size + (1:local_mat_size);
        range_test = mesh_test.T((i_vec_test-1)*mesh_test.N_lb+1:i_vec_test*mesh_test.N_lb,n);
        range_trial = mesh_trial.T((i_vec_trial-1)*mesh_trial.N_lb+1:i_vec_trial*mesh_trial.N_lb,n)';
        i = repmat(range_test,mesh_trial.N_lb,1);
        j = repmat(range_trial,mesh_test.N_lb,1);
        j = j(:);
        index_i(count) = i;
        index_j(count) = j;
        index_v(count) = mat(:);

    end
end

A = sparse(index_i,index_j,index_v,matrix_size(1),matrix_size(2));

