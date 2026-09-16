function solver_2D_Poisson_d
clear;

Gauss_type_2D = 3;

tic;
[T,P] = mesh_square_triangle(-1,1,-1,1,100,100);
mesh = mesh_topo(T,P);
hmax = mesh.hmax
mesh_trial = mesh_FE('P1',1,mesh,Gauss_type_2D,0,[]);
mesh_test = mesh_trial;

coe_fun = @fun_c;
A = assemble_matrix_2D(coe_fun,[1,1],mesh,mesh_trial,mesh_test,[1,1,1,0,1,0;1,1,0,1,0,1],Gauss_type_2D);

rhs_fun = @fun_rhs;
b = assemble_vector_2D(rhs_fun,mesh,mesh_test,0,0,1,Gauss_type_2D);

bndy_fun = @fun_bndy;
[A,b] = treat_Dirichlet(A,b,mesh,mesh_trial,0,bndy_fun);

solution = A\b;
toc;

exact_fun = @fun_exact;
print_all_error(solution,exact_fun,mesh,mesh_trial,Gauss_type_2D);
