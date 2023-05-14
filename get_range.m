function result = get_range(mesh_FE,i_elem,i_vec)
    % get dof indices of a local element.
    
range = (i_vec-1)*mesh_FE.N_lb+(1:mesh_FE.N_lb);
result = mesh_FE.T(range,i_elem);