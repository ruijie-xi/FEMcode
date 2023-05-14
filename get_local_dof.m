function local_dof = get_local_dof(u,~,mesh_FE,i_elem,i_vec)
% get dof indices of a local element.

range = (i_vec-1)*mesh_FE.N_lb+(1:mesh_FE.N_lb);
local_dof = u(mesh_FE.T(range,i_elem));