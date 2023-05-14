function result = basis_local_read(mesh_FE,i_elem,i_vec,dx,dy)
% read local basis information from cached data

range = (i_elem-1)*mesh_FE.N_lb + (1:mesh_FE.N_lb);
if dx==0 && dy==0
    result = mesh_FE.basis(i_vec,range,:);
elseif dx==1 && dy==0
    result = mesh_FE.basis_dx(i_vec,range,:);
elseif dx==0 && dy==1
    result = mesh_FE.basis_dy(i_vec,range,:);
end

result = reshape(result,mesh_FE.N_lb,[]);