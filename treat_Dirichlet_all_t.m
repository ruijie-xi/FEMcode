function [A,b] = treat_Dirichlet_all_t(A,b,mesh,mesh_FE,Dbndy_info,t)


for i_D = 1:size(Dbndy_info,1)
    % {bndy_number, bndy_type, bndy_info(bndy_fun or dof_on_edge)}
    if Dbndy_info{i_D,2} == 0 % Dirichlet boundary
        bndy_fun = Dbndy_info{i_D,3};
        [A,b] = treat_Dirichlet(A,b,mesh,mesh_FE,Dbndy_info{i_D,1},@(x)bndy_fun(x,t));
    elseif Dbndy_info{i_D,2} == 1 % set some dofs to zero. (like u\cdot n = 0)
        [A,b] = treat_Dirichlet_zero(A,b,mesh,mesh_FE,Dbndy_info{i_D,1},Dbndy_info{i_D,3});
    end
end

end