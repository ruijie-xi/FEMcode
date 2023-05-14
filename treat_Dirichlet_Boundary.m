function [A,b] = treat_Dirichlet_Boundary(A0,b0,mesh,mesh_FE,bndy_fun)
% apply Dirichlet Boundary conditions on linear system.

[M,N] = size(A0);
large = 1e8;

N_bndy0 = size(mesh_FE.Dbndynodes,2);

index_b_i = zeros(N_bndy0,1);
index_b_j = ones(N_bndy0,1);
index_b_v = zeros(N_bndy0,1);

index_A_i = zeros(N_bndy0,1);
index_A_j = zeros(N_bndy0,1);
index_A_v = zeros(N_bndy0,1);

for i = 1:N_bndy0
    idx = mesh_FE.Dbndynodes(i);
    index_b_i(i) = idx;
    index_b_v(i) = large*compute_dof(bndy_fun,mesh,mesh_FE,idx)-b0(idx,1);
    index_A_i(i) = idx;
    index_A_j(i) = idx;
    index_A_v(i) = large-A0(idx,idx);
end

Ap = sparse(index_A_i,index_A_j,index_A_v,M,N);
bp = sparse(index_b_i,index_b_j,index_b_v,M,1);

A = A0+Ap;
b = b0+bp;