function err = compute_norm(vector,mesh,mesh_trial,error_type,Gauss_type)
% compute norms of FE functions.

if strcmp(error_type,'L2')
    err = 0;
    for n = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,n));
        [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        for i_vec = 1:mesh_trial.dim
            err = err + Gauss_quad_2D_norm(vector,pts,weights,...
                mesh,mesh_trial,n,0,0,i_vec);
        end
    end
    err = sqrt(err);
elseif strcmp(error_type,'H1')
    err = 0;
    for n = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,n));
        [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        for i_vec = 1:mesh_trial.dim
            err = err + Gauss_quad_2D_norm(vector,pts,weights,...
                mesh,mesh_trial,n,0,0,i_vec);
            err = err + Gauss_quad_2D_norm(vector,pts,weights,...
                mesh,mesh_trial,n,1,0,i_vec);
            err = err + Gauss_quad_2D_norm(vector,pts,weights,...
                mesh,mesh_trial,n,0,1,i_vec);
        end
    end
    err = sqrt(err);

elseif strcmp(error_type,'H10')
    err = 0;
    for n = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,n));
        [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        for i_vec = 1:mesh_trial.dim
            err = err + Gauss_quad_2D_norm(vector,pts,weights,...
                mesh,mesh_trial,n,1,0,i_vec);
            err = err + Gauss_quad_2D_norm(vector,pts,weights,...
                mesh,mesh_trial,n,0,1,i_vec);
        end
    end
    err = sqrt(err);
elseif strcmp(error_type,'L_inf')
    err = 0;
    for n = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,n));
        [~,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        for i_vec = 1:mesh_trial.dim
            err = max([err,Gauss_Linf_norm(vector,pts,mesh,mesh_trial,n,i_vec)]);
        end
    end
end


