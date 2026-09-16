function err = compute_error(solution,u_fun,mesh,mesh_trial,error_type,Gauss_type)
% compute error |solution - u_fun|.
% error_type:
%   'L2','H1','H10'(seminorm),'L_inf','Hdiv'.
%   'Hdiv' is the full norm: sqrt(||e||_L2^2 + ||div(e)||_L2^2).

if strcmp(error_type,'L2')
    err = 0;
    for n = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,n));
        [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        for i_vec = 1:mesh_trial.dim
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
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
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
                mesh,mesh_trial,n,0,0,i_vec);
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
                mesh,mesh_trial,n,1,0,i_vec);
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
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
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
                mesh,mesh_trial,n,1,0,i_vec);
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
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
            err = max([err,Gauss_Linf_error(u_fun,solution,pts,mesh,mesh_trial,n,i_vec)]);
        end
    end
elseif strcmp(error_type,'Hdiv')
    assert(mesh_trial.dim==2);
    err = 0;
    for n = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,n));
        [weights,pts] = generate_Gauss_local_triangle(vertices,Gauss_type);
        for i_vec = 1:mesh_trial.dim
            err = err + Gauss_quad_2D_error(u_fun,solution,pts,weights,...
                mesh,mesh_trial,n,0,0,i_vec);
        end
        exact_dx = u_fun(pts,1,0);
        exact_dy = u_fun(pts,0,1);
        div_error = exact_dx(1,:) + exact_dy(2,:) ...
                  - FE_function_local_2D_read(solution,mesh,mesh_trial,n,1,1,0) ...
                  - FE_function_local_2D_read(solution,mesh,mesh_trial,n,2,0,1);
        err = err + div_error.^2*weights';
    end
    err = sqrt(err);
end


