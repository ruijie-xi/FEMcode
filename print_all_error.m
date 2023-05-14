function error = print_all_error(solution,u_fun,mesh,mesh_trial,Gauss_type)

L2error = compute_error(solution,u_fun,mesh,mesh_trial,'L2',Gauss_type);
H1error = compute_error(solution,u_fun,mesh,mesh_trial,'H1',Gauss_type);
Linferror = compute_error(solution,u_fun,mesh,mesh_trial,'L_inf',Gauss_type);
fprintf('L2 error: %.2e\n',L2error);
fprintf('H1 error: %.2e\n',H1error);
fprintf('L_inf error: %.2e\n',Linferror);

error.L2 = L2error;
error.H1 = H1error;
error.Linf = Linferror;