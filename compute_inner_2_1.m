function val = compute_inner_2_1(mesh,vec_1_1,mesh_FE_1_1,vec_1_2,mesh_FE_1_2,vec_2,mesh_FE_2,inner_info,Gauss_type)
% compute inner product.

val = 0;
num = size(inner_info,1);
for n = 1:num
    i_vec_1_1 = inner_info(n,1);
    dx_1_1 = inner_info(n,2);
    dy_1_1 = inner_info(n,3);
    i_vec_1_2 = inner_info(n,4);
    dx_1_2 = inner_info(n,5);
    dy_1_2 = inner_info(n,6);
    i_vec_2 = inner_info(n,7);
    dx_2 = inner_info(n,8);
    dy_2 = inner_info(n,9);
    for i_elem = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,i_elem));
        [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        val = val + Gauss_quad_2D_inner_2_1(mesh,i_elem,vec_1_1,mesh_FE_1_1,dx_1_1,dy_1_1,i_vec_1_1,...
            vec_1_2,mesh_FE_1_2,dx_1_2,dy_1_2,i_vec_1_2,vec_2,mesh_FE_2,dx_2,dy_2,i_vec_2,pts,weights);
    end
end