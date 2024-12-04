function [A,b] = treat_Dirichlet(A0,b0,mesh,mesh_FE,Dbndy_number,bndy_fun)

[M,N] = size(A0);
large = 1e0;

Dbndynodes = get_Dbndynodes(mesh,mesh_FE,Dbndy_number);
N_bndy0 = numel(Dbndynodes);

A0(Dbndynodes,:) = 0;
b0(Dbndynodes) = 0;

index_b_i = zeros(N_bndy0,1);
index_b_j = ones(N_bndy0,1);
index_b_v = zeros(N_bndy0,1);

index_A_i = zeros(N_bndy0,1);
index_A_j = zeros(N_bndy0,1);
index_A_v = zeros(N_bndy0,1);

for i = 1:N_bndy0
    idx = Dbndynodes(i);
    index_b_i(i) = idx;
    index_b_v(i) = large*compute_dof(bndy_fun,mesh,mesh_FE,idx)-b0(idx,1);
    
    index_A_i(i) = idx;
    index_A_j(i) = idx;
    index_A_v(i) = large;

end

Ap = sparse(index_A_i,index_A_j,index_A_v,M,N);
bp = sparse(index_b_i,index_b_j,index_b_v,M,1);

A = A0+Ap;
b = b0+bp;