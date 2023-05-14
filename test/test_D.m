Gauss_type_2D = 3;

[T,P] = mesh_square_triangle(-1,1,-1,1,100,100);
mesh = mesh_topo(T,P);
hmax = mesh.hmax
mesh_trial = mesh_FE('Ned2-1',2,mesh,Gauss_type_2D,0,[]);
mesh_test = mesh_trial;

bndy_fun = @(x)0;
[A,b] = treat_Dirichlet(0,0,mesh,mesh_trial,0,bndy_fun);
