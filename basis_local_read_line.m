function result = basis_local_read_line(mesh,mesh_FE,i_elem,i_edge,i_vec,dx,dy)
% read local basis value (on edges) from cached data.

% i = find(abs(mesh.TE(:,i_elem))==i_edge,1);

i = mesh.ET(2,i_edge);
if(abs(mesh.TE(i,i_elem))~=i_edge)
    i = mesh.ET(1,i_edge);
end

range = (i_elem-1)*3*mesh_FE.N_lb + (i-1)*mesh_FE.N_lb + (1:mesh_FE.N_lb);
if dx==0 && dy==0
    result = mesh_FE.basis_line(i_vec,range,:);
elseif dx==1 && dy==0
    result = mesh_FE.basis_line_dx(i_vec,range,:);
elseif dx==0 && dy==1
    result = mesh_FE.basis_line_dy(i_vec,range,:);
end

result = reshape(result,mesh_FE.N_lb,[]);