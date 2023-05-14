function val = compute_integral_FE(mesh,vec_FE,mesh_FE,dx,dy,i_vec,Gauss_type)
%compute integral of an FE function.

val = 0;
for n = 1:mesh.N_elem
    vertices = mesh.P(:,mesh.T(:,n));
    [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
    val = val + Gauss_quad_2D_FE(vec_FE,pts,weights,mesh,mesh_FE,n,dx,dy,i_vec);
end
