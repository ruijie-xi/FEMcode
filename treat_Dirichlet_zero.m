function [A,b] = treat_Dirichlet_zero(A0,b0,mesh,mesh_FE,Dbndy_number,Dof_on_edge)

[M,N] = size(A0);
large = 1e8;

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
N_bndy0 = numel(Dbndynodes);

index_b_i = zeros(N_bndy0,1);
index_b_j = ones(N_bndy0,1);
index_b_v = zeros(N_bndy0,1);

index_A_i = zeros(N_bndy0*N,1);
index_A_j = zeros(N_bndy0*N,1);
index_A_v = zeros(N_bndy0*N,1);

for i = 1:N_bndy0
    idx = Dbndynodes(i);
    index_b_i(i) = idx;
    index_b_v(i) = -b0(idx,1);
    for j = 1:N
        index_A_i((i-1)*N+j) = idx;
        index_A_j((i-1)*N+j) = j;
        if j==idx
            index_A_v((i-1)*N+j) = large-A0(idx,idx);
        else
            index_A_v((i-1)*N+j) = -A0(idx,j);
        end
    end
end


Ap = sparse(index_A_i,index_A_j,index_A_v,M,N);
bp = sparse(index_b_i,index_b_j,index_b_v,M,1);

A = A0+Ap;
b = b0+bp;