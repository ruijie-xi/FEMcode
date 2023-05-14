function vec = local_vector_line(i_edge,coe_fun,mesh,mesh_test,...
    i_vec_test,dx_test,dy_test,Gauss_type)
% assemble local vector (integral on line)

i_elem = mesh.E(2,i_edge);
endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
[Gauss_weights,Gauss_pts] = generate_Gauss_local_line(endpts,Gauss_type);

vec = basis_local_read_line(mesh,mesh_test,i_elem,i_edge,i_vec_test,dx_test,dy_test)*...
    (Gauss_weights.*coe_fun(Gauss_pts))';

