function A = assemble_matrix_2D_FE_2(coe_fun,coe_num,mesh,mesh_trial,mesh_test,...
    vec_FE_1,mesh_FE_1,vec_FE_2,mesh_FE_2,mat_info,Gauss_type)
% assemble matrix with 2 FE functions.

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
        i_vec_FE_1 = mat_info(i_mat,7);
        dx_FE_1 = mat_info(i_mat,8);
        dy_FE_1 = mat_info(i_mat,9);
        i_vec_FE_2 = mat_info(i_mat,10);
        dx_FE_2 = mat_info(i_mat,11);
        dy_FE_2 = mat_info(i_mat,12);
        
        mat = coe_num(i_mat)*local_matrix_FE_FE(n,coe_fun,mesh,mesh_trial,mesh_test,i_vec_trial,dx_trial,dy_trial,i_vec_test,dx_test,dy_test,vec_FE_1,mesh_FE_1,i_vec_FE_1,dx_FE_1,dy_FE_1,vec_FE_2,mesh_FE_2,i_vec_FE_2,dx_FE_2,dy_FE_2,Gauss_type);
    
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

