function Dbndynodes = get_Dbndynodes(mesh,mesh_FE,Dbndy_number)

tmp = 1:mesh.N_edge;
flag = mesh.E(1,:)==Dbndy_number;
Dbndy = tmp(flag);
N_Dbndy = numel(Dbndy);
Dbndynodes = zeros(1,N_Dbndy*size(mesh_FE.Edgedofs,2));

for i_bndyedge = 1:numel(Dbndy)
    i_edge = Dbndy(i_bndyedge);
    i_elem = mesh.E(2,i_edge);
    i_edge_elem = mesh.ET(2,i_edge);
    Edgedofs = mesh_FE.Edgedofs(i_edge_elem,:);
    Globaldofs = mesh_FE.T(Edgedofs,i_elem);
    Dbndynodes((i_bndyedge-1)*size(mesh_FE.Edgedofs,2)+1:i_bndyedge*size(mesh_FE.Edgedofs,2)) = Globaldofs;
end

Dbndynodes = unique(Dbndynodes);