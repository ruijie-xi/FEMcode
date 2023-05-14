function plot_vec(vector,mesh,mesh_FE,varargin)
% plot a vector-valued function.

if nargin==3
    dx = 0;
    dy = 0;
else
    dx = varargin{1};
    dy = varargin{2};
end

n_plot = 4;
X = zeros(1,mesh.N_elem*n_plot);
Y = zeros(1,mesh.N_elem*n_plot);
U = zeros(1,mesh.N_elem*n_plot);
V = zeros(1,mesh.N_elem*n_plot);

for i_elem = 1:mesh.N_elem
    vertices = mesh.P(:,mesh.T(:,i_elem));
    drawing_pts = [vertices,mean(vertices,2)];

    X((i_elem-1)*n_plot+(1:n_plot)) = drawing_pts(1,:);
    Y((i_elem-1)*n_plot+(1:n_plot)) = drawing_pts(2,:);
    
    U((i_elem-1)*n_plot+(1:n_plot)) = ...
        FE_function_local_2D(drawing_pts,vector,mesh,mesh_FE,i_elem,1,dx,dy);
    V((i_elem-1)*n_plot+(1:n_plot)) = ...
        FE_function_local_2D(drawing_pts,vector,mesh,mesh_FE,i_elem,2,dx,dy);
end

quiver(X,Y,U,V);