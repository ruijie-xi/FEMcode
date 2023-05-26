function u = set_bndy_zero(u,mesh,mesh_FE,Dbndy_number,Dof_on_edge)
% set some dofs on bndy to zero.

tmp = 1:mesh.N_edge;
flag = mesh.E(1,:)==Dbndy_number;
Dbndy = tmp(flag);
N_Dbndy = numel(Dbndy);
Dbndynodes = zeros(1,N_Dbndy*numel(Dof_on_edge));

for i_bndyedge = 1:numel(Dbndy)
    i_edge = Dbndy(i_bndyedge);
    i_elem = mesh.E(2,i_edge);
    i_edge_elem = mesh.ET(2,i_edge);
    Edgedofs = mesh_FE.Edgedofs(i_edge_elem,Dof_on_edge);
    Globaldofs = mesh_FE.T(Edgedofs,i_elem);
    Dbndynodes((i_bndyedge-1)*numel(Dof_on_edge)+1:i_bndyedge*numel(Dof_on_edge)) = Globaldofs;
end

Dbndynodes = unique(Dbndynodes);

u(Dbndynodes) = 0;